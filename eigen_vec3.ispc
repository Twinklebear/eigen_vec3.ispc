#define DIM 3
#define MAX_ITER 30
#define SQRT3 1.73205080756887729352744634151f
#define INV3 0.33333333333333333333333333333f
// Access row r and column c of a row-major 3x3 matrix
#define AT(m, r, c) ((m)[(r) * 3 + (c)])
#define EPSILON 1.19209290E-07f
// Select between analytic only or analytic + QL fallback if potentially large error may occur
#define ANALYTIC_ONLY 1

#ifndef sqr
// Will this be ok with integrating with ospray's built in sqr function?
inline float sqr(float x){
	return x * x;
}
#endif

// Calculates the eigenvalues of a symmetric 3x3 matrix `mat` using Cardano's
// analytical algorithm. The matrix should contain all entries of the matrix for now
// but TODO we only really need the upper triangle.
// `w` will contain the computed eigenvalues
// This is a port of dsyevc3 from the Hybrid solver described by Joachim Kopp
// See: http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
// TODO: can we do something about the required scatter/gathers to load and store the results?
void eig_vals(const uniform float *m, float *w){
	// Determine coefficients of characteristic poynomial. We write
	//       | a   d   f  |
	//  A =  | d*  b   e  |
	//       | f*  e*  c  |
	float de = AT(m, 0, 1) * AT(m, 1, 2);
	float dd = sqr(AT(m, 0, 1));
	float ee = sqr(AT(m, 1, 2));
	float ff = sqr(AT(m, 0, 2));
	float trace = AT(m, 0, 0) + AT(m, 1, 1) + AT(m, 2, 2);
	float c1 = (AT(m, 0, 0) * AT(m, 1, 1) + AT(m, 0, 0) * AT(m, 2, 2) + AT(m, 1, 1) * AT(m, 2, 2))
		- (dd + ee + ff);
	float c0 = AT(m, 2, 2) * dd + AT(m, 0, 0) * ee + AT(m, 1, 1) * ff
		- AT(m, 0, 0) * AT(m, 1, 1) * AT(m, 2, 2) - 2.f * AT(m, 0, 2) * de;

	float p = sqr(trace) - 3.f * c1;
	float q = trace * (p - (3.f / 2.f) * c1) - (27.f / 2.f) * c0;
	float sqrt_p = sqrt(abs(p));

	float phi = 27.f * (0.25f * sqr(c1) * (p - c1) + c0 * (q + 27.f / 4.f * c0));
	phi = (1.f / 3.f) * atan2(sqrt(abs(phi)), q);

	float c = sqrt_p * cos(phi);
	float s = (1.f / SQRT3) * sqrt_p * sin(phi);

	w[1] = (1.f / 3.f) * (trace - c);
	w[2] = w[1] + s;
	w[0] = w[1] + c;
	w[1] -= s;
}
// Reduce a symmetric 3x3 matrix to tridiagonal form by applying (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// A is the 3x3 matrix to be diagonalized
// Q will be used to store the eigenvectors and should also be 3x3
// d will store the diagonal and should have 3 elements
// e will store the off-diagonal and should have 2 elements
// This is a port of dsytrd3 from the Hybrid solver described by Joachim Kopp
// See: http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
void reduce_tridiag(const uniform float *A, float *Q, float *d, varying float * varying e){
	// Initialize Q to the identity matrix
	memset(Q, 0, sizeof(float) * 9);
	AT(Q, 0, 0) = 1.f;
	AT(Q, 1, 1) = 1.f;
	AT(Q, 2, 2) = 1.f;

	// Bring first row and col to desired form
	float h = sqr(AT(A, 0, 1)) + sqr(AT(A, 0, 2));
	float g = sqrt(h);
	if (AT(A, 0, 1) > 0.f){
		g = -g;
	}
	e[0] = g;
	float f = g * AT(A, 0, 1);
	float omega = h - f;

	float u[3] = { 0, AT(A, 0, 1) - g, AT(A, 0, 2) };
	float q[3] = { 0 };
	if (omega > 0.f){
		omega = 1.f / omega;
		float kappa = 0.f;
		for (int i = 1; i < DIM; ++i){
			f = AT(A, 1, i) * u[1] + AT(A, i, 2) * u[2];
			q[i] = omega * f;
			kappa += u[i] * f;
		}

		kappa *= 0.5f * sqr(omega);
		for (int i = 1; i < DIM; ++i){
			q[i] = q[i] - kappa * u[i];
		}
		d[0] = AT(A, 0, 0);
		d[1] = AT(A, 1, 1) - 2.f * q[1] * u[1];
		d[2] = AT(A, 2, 2) - 2.f * q[2] * u[2];

		// Store inverse Householder transform in Q
		for (int j = 1; j < DIM; ++j){
			f = omega * u[j];
			for (int i = 1; i < DIM; ++i){
				AT(Q, i, j) = AT(Q, i, j) - f * u[i];
			}
		}
		// Calculate update A[1][2] and store in e[1]
		e[1] = AT(A, 1, 2) - q[1] * u[2] - u[1] * q[2];
	}
	else {
		for (int i = 0; i < DIM; ++i){
			d[i] = AT(A, i, i);
		}
		e[1] = AT(A, 1, 2);
	}
}
// Calculates eigenvalues and normalized eigenvectors of a symmetric 3x3 matrix A
// using the QL algorithm with implicit shifts, preceded by a Householder reduction
// to tridiagonal form.
// Q should be a 3x3 matrix and will store the 3 eigenvectors
// w will store the 3 eigenvalues
// returns 0 on success, -1 if we didn't converge to a solution
// This is a port of dsyevq3 from the Hybrid solver described by Joachim Kopp
// See: http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
int eig_ql(const uniform float *A, float *Q, float *w){
	// Third element used as temp storage
	float e[3] = { 0 };
	// Transform A to a real tridiagonal form using Householder method
	reduce_tridiag(A, Q, w, &e[0]);

	int m = 0;
	float g, r, p, f, b, s, c, t;
	// Calculate eigensystem of remaining real symmetric tridiagonal matrix
	// with the QL method
	for (int l = 0; l < DIM - 1; ++l){
		int n_iter = 0;
		while (1){
			// Check for convergence and exit if off-diagonal elem e(l) is 0
			for (m = l; m <= DIM - 2; ++m){
				g = abs(w[m]) + abs(w[m + 1]);
				if (abs(e[m]) + g == g){
					break;
				}
			}
			if (m == l){
				break;
			}
			if (++n_iter > MAX_ITER){
				return -1;
			}

			g = (w[l + 1] - w[l]) / (e[l] + e[l]);
			r = sqrt(sqr(g) + 1.f);
			if (g > 0.f){
				g = w[m] - w[l] + e[l] / (g + r);
			}
			else {
				g = w[m] - w[l] + e[l] / (g - r);
			}

			s = 1.f;
			c = 1.f;
			p = 0.f;
			for (int i = m - 1; i >= l; --i){
				f = s * e[i];
				b = c * e[i];
				if (abs(f) > abs(g)){
					c = g / f;
					r = sqrt(sqr(c) + 1.f);
					e[i + 1] = f * r;
					s = rcp(r);
					c *= s;
				}
				else {
					s = f / g;
					r = sqrt(sqr(s) + 1.f);
					e[i + 1] = g * r;
					c = rcp(r);
					s *= c;
				}

				g = w[i + 1] - p;
				r = (w[i] - g) * s + 2.f * c * b;
				p = s * r;
				w[i + 1] = g + p;
				g = c * r - b;

				// Form eigenvectors
				for (int k = 0; k < DIM; ++k){
					t = AT(Q, k, i + 1);
					AT(Q, k, i + 1) = s * AT(Q, k, i) + c * t;
					AT(Q, k, i) = c * AT(Q, k, i) - s * t;
				}
			}
			w[l] -= p;
			e[l] = g;
			e[m] = 0.f;
		}
	}
	return 0;
}
// Calculates the eigenvalues and normalized eigenvectors of the symmetric 3x3
// matrix A using Cardano's method for eigenvalue and analytical method for
// the eigenvectors. If it's determined that a large amount of error will
// occur in the results using this method it falls back to the iterative QL algorithm
// TODO: We only use the diagonal and upper-triangle of A so we really don't need the full 9 floats
// Q: contains the eigenvectors in the columns, is stored row-major
// w: the eigenvalues
// returns -1 if the system couldn't be solved
// This is a port of dsyevh3 from the Hybrid solver described by Joachim Kopp
// See: http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/
int eig_hybrid(const uniform float *A, float *Q, float *w){
	// Compute the eigenvalues
	eig_vals(A, w);

	float t = abs(w[0]);
	float u = abs(w[1]);
	if (u > t){
		t = u;
	}
	u = abs(w[2]);
	if (u > t){
		t = u;
	}
	if (t < 1.f){
		u = t;
	}
	else {
		u = sqrt(t);
	}

	AT(Q, 0, 1) = AT(A, 0, 1) * AT(A, 1, 2) - AT(A, 0, 2) * AT(A, 1, 1);
	AT(Q, 1, 1) = AT(A, 0, 2) * AT(A, 0, 1) - AT(A, 1, 2) * AT(A, 0, 0);
	AT(Q, 2, 1) = sqr(AT(A, 0, 1));

	// Compute first eigenvector using
	// v[0] = (A - w[0]).eq x (A - w[0]).e2
	AT(Q, 0, 0) = AT(Q, 0, 1) + AT(A, 0, 2) * w[0];
	AT(Q, 1, 0) = AT(Q, 1, 1) + AT(A, 1, 2) * w[0];
	AT(Q, 2, 0) = (AT(A, 0, 0) - w[0]) * (AT(A, 1, 1) - w[0]) - AT(Q, 2, 1);
	float norm = sqr(AT(Q, 0, 0)) + sqr(AT(Q, 1, 0)) + sqr(AT(Q, 2, 0));

	// If vectors are nearly linearly dependent or if there may have been large
	// cancellation in A[i][i] - w[0] fall back to QL algorithm.
#ifndef ANALYTIC_ONLY
	const float error = 256.f * EPSILON * sqr(u);
	if (norm <= error){
		return eig_ql(A, Q, w);
	}
	else {
#endif
		norm = sqrt(rcp(norm));
		for (int j = 0; j < DIM; ++j){
			AT(Q, j, 0) = AT(Q, j, 0) * norm;
		}
#ifndef ANALYTIC_ONLY
	}
#endif

	// Calculate second eigenvector using
	// v[1] = (A - w[1]).e1 x (A - w[1]).e2
	AT(Q, 0, 1) = AT(Q, 0, 1) + AT(A, 0, 2) * w[1];
	AT(Q, 1, 1) = AT(Q, 1, 1) + AT(A, 1, 2) * w[1];
	AT(Q, 2, 1) = (AT(A, 0, 0) - w[1]) * (AT(A, 1, 1) - w[1]) - AT(Q, 2, 1);
	norm = sqr(AT(Q, 0, 1)) + sqr(AT(Q, 1, 1)) + sqr(AT(Q, 2, 1));
#ifndef ANALYTIC_ONLY
	if (norm <= error){
		return eig_ql(A, Q, w);
	}
	else {
#endif
		norm = sqrt(rcp(norm));
		for (int j = 0; j < DIM; ++j){
			AT(Q, j, 1) = AT(Q, j, 1) * norm;
		}
#ifndef ANALYTIC_ONLY
	}
#endif

	// Calculate the third eigenvector using cross product
	AT(Q, 0, 2) = AT(Q, 1, 0) * AT(Q, 2, 1) - AT(Q, 2, 0) * AT(Q, 1, 1);
	AT(Q, 1, 2) = AT(Q, 2, 0) * AT(Q, 0, 1) - AT(Q, 0, 0) * AT(Q, 2, 1);
	AT(Q, 2, 2) = AT(Q, 0, 0) * AT(Q, 1, 1) - AT(Q, 1, 0) * AT(Q, 0, 1);
	return 0;
}
// Calculates eigenvectors and eigenvalues of a symmetric 3x3 matrix using
// the analytic method from Jeong '06
// See: http://people.seas.harvard.edu/~wkjeong/publication/vg06-jeong.pdf
// TODO: Is this implemented properly? It seems to give much less accurate eigenvectors
// it also isn't that much faster so maybe just go with the hybrid method?
int eig_jeong(const uniform float *A, float *Q, float *w){
	// Compute eigenvalues
	eig_vals(A, w);
	// Compute eigenvectors
	float<3> va = { AT(A, 0, 0) - w[0], AT(A, 0, 0) - w[1], AT(A, 0, 0) - w[2] };
	float<3> vb = { AT(A, 1, 1) - w[0], AT(A, 1, 1) - w[1], AT(A, 1, 1) - w[2] };
	float<3> vc = { AT(A, 2, 2) - w[0], AT(A, 2, 2) - w[1], AT(A, 2, 2) - w[2] };

	float p = AT(A, 0, 1) * AT(A, 1, 2);
	float<3> vd = { p - vb.x * AT(A, 0, 2), p - vb.y * AT(A, 0, 2), p - vb.z * AT(A, 0, 2) };

	p = AT(A, 0, 2) * AT(A, 1, 2);
	float<3> ve = { p - vc.x * AT(A, 0, 1), p - vc.y * AT(A, 0, 1), p - vc.z * AT(A, 0, 1) };

	p = AT(A, 0, 2) * AT(A, 0, 1);
	float<3> vf = { p - va.x * AT(A, 1, 2), p - va.y * AT(A, 1, 2), p - va.z * AT(A, 1, 2) };

	va = vd * ve;
	vb = ve * vf;
	vc = vd * vf;

	float<3> ev1 = { va.x, vb.x, vc.x };
	float<3> ev2 = { va.y, vb.y, vc.y };
	float<3> ev3 = { va.z, vb.z, vc.z };
	float norm = sqrt(rcp(sqr(ev1.x) + sqr(ev1.y) + sqr(ev1.z)));
	ev1 = ev1 * norm;

	norm = sqrt(rcp(sqr(ev2.x) + sqr(ev2.y) + sqr(ev2.z)));
	ev2 = ev2 * norm;

	norm = sqrt(rcp(sqr(ev3.x) + sqr(ev3.y) + sqr(ev3.z)));
	ev3 = ev3 * norm;

	AT(Q, 0, 0) = ev1.x;
	AT(Q, 1, 0) = ev1.y;
	AT(Q, 2, 0) = ev1.z;

	AT(Q, 0, 1) = ev2.x;
	AT(Q, 1, 1) = ev2.y;
	AT(Q, 2, 1) = ev2.z;

	AT(Q, 0, 2) = ev3.x;
	AT(Q, 1, 2) = ev3.y;
	AT(Q, 2, 2) = ev3.z;
	return 0;
}

export void eig_vec_test(const uniform float matrices[], uniform float vecs[], uniform float vals[],
        const uniform int num_mats){
    foreach (i = 0 ... num_mats){
        eig_hybrid(&matrices[i * 9], &vecs[i * 9], &vals[i * 3]);
    }
}

