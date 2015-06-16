/*
 * Uses a Roe solver to solve the Burger equation. The graphing depends on StdDraw.java, which may be found here.
 */

public class Burgers {

	double[] q;
	//toggles Gaussian or shock tube initial data
	boolean smooth;
	//number of time steps
	int numsteps;
	//length of time steps
	double delta_t;

	//number of cells
	int Nx;
	//number of ghost cells at each boundary
	int Ng;

	//domain
	double xmin;
	double xmax;
	//distance between cells
	double delta_x;
	//holds positions of cell centers
	double[] x;


	public Burgers(int Nx, double dt, boolean smooth) {
		this.smooth = smooth;
		delta_t = dt;
		this.Nx = Nx;

		xmin = 0;
		xmax = 5;
		Ng = 2;

		delta_x = (xmax - xmin)/Nx;
		x = new double[Nx];
		q = new double[Nx];
		for (int i = 0; i < Nx; i++) {
			x[i] = xmin + (i + .5)*delta_x;
		}

	}

	public Burgers(double xmin, double xmax, int Nx, double dt, boolean smooth) {
		this.smooth = smooth;
		delta_t = dt;
		this.Nx = Nx;

		this.xmin = xmin;
		this.xmax = xmax;
		Ng = 2;

		delta_x = (xmax - xmin)/Nx;
		x = new double[Nx];
		q = new double[Nx];
		for (int i = 0; i < Nx; i++) {
			x[i] = xmin + (i + .5)*delta_x;
		}


	}

	public void init_fluid() {

		// if smooth = true, Gaussian initial data, otherwise shock tube
		if (smooth == true) {
			for (int i = 0; i < Nx; i++) {
				q[i] = Math.exp(-(x[i] - (xmax + xmin)/2)*(x[i] - (xmax + xmin)/2));
			}
		}

		else {
			//shock tube
			double q_L = 0.1;
			double q_R = 1;
			for (int i = 0; i < Nx; i++) {
				if (i < Nx/2) q[i] = q_L;
				else q[i] = q_R;
			}
		}
	}

	//minmod function
	private double minmod(double a, double b){
		if (a*b < 0) return 0;
		else if(Math.abs(a) < Math.abs(b)) return a;
		else return b;
	}

	//calculating Jacobian matrix at x
	private double calc_A(double q) {
		//the Jacobian is just q
		double A = q;
		return A;
	}

	//calculating eigenvectors of Jacobian
	private double calc_eta(double A) {
		//scalar, so eigenvector is just 1
		double eta = 1;
		return eta;
	}

	//calculating eigenvalues of Jacobian
	private double calc_lambda(double A) {
		//just the same as A in the scalar case
		double lambda = A;
		return lambda;
	}

	//calculating jumps
	private double calc_omega(double qL, double qR, double eta) {
		double omega = qR - qL;
		return omega;
	}

	//calculate Roe flux
	private double calc_FF(double fR, double fL, double lambda, double eta, double omega) {
		double FF = 0.5*(fR + fL - lambda*eta*omega);
		return FF;
	}

	//flux as a function of q
	private double calc_f(double q) {
		double f = 0.5*q*q;
		return f;
	}
	
	private double dqL_1(double dq_imh, double dq_iph) {
		return minmod(dq_iph, dq_imh);
	}
	
	private double dqR_1(double dq_iph, double dq_ip3h) {
		return minmod(dq_iph, dq_ip3h);
	}
	
	private double dqL_2(double dq_imh, double dq_iph) {
		if (Math.abs(dq_iph) <= 1e-6) return 0;
		else {
			double theta = dq_imh/dq_iph;
			return Math.max(Math.max(0, Math.min(1, 2*theta)), Math.min(2, theta))*dq_iph;
		}
	}
	
	private double dqR_2(double dq_iph, double dq_ip3h) {
		if (Math.abs(dq_iph) <= 1e-6) return 0;
		else {
			double theta = dq_ip3h/dq_iph;
			return Math.max(Math.max(0, Math.min(1, 2*theta)), Math.min(2, theta))*dq_iph;
		}
	}
	
	private double dqL_3(double dq_imh, double dq_iph) {
		if (Math.abs(dq_iph) <= 1e-6) return 0;
		else {
			double theta = dq_imh/dq_iph;
			return Math.max(0, Math.min(Math.min((1 + theta)/2, 2), 2*theta))*dq_iph;
		}
	}
	
	private double dqR_3(double dq_iph, double dq_ip3h) {
		if (Math.abs(dq_iph) <= 1e-6) return 0;
		else {
			double theta = dq_ip3h/dq_iph;
			return Math.max(0, Math.min(Math.min((1 + theta)/2, 2), 2*theta))*dq_iph;
		}
	}
	
	private double dqL_4(double s_imh, double s_iph) {
		if (Math.abs(s_iph) <= 1e-6) return 0;
		else {
			double theta = s_imh/s_iph;
			return s_iph*(theta + Math.abs(theta))/(1 + Math.abs(theta));
		}
	}
	
	private double dqR_4(double s_iph, double s_ip3h) {
		if (Math.abs(s_iph) <= 1e-6) return 0;
		else {
			double theta = s_ip3h/s_iph;
			return s_iph*(theta + Math.abs(theta))/(1 + Math.abs(theta));
		}
	}
	
	//reconstructs primitive variable q from the left
	private double recos_qL(double q_im1, double q_i, double q_ip1, double r_im1, double r_i, double r_ip1) {
		double dq_iph = q_ip1 - q_i;
		double dq_imh = q_i - q_im1;
		double dq_i = dqL_2(dq_imh, dq_iph);
		double rqL_iph = q_i + dq_i*0.5;
		return rqL_iph;
	}

	//reconstructs primitive variable q from the right
	private double recos_qR(double q_i, double q_ip1, double q_ip2, double r_i, double r_ip1, double r_ip2) {
		double dq_iph = q_ip1 - q_i;
		double dq_ip3h = q_ip2 - q_ip1;
		double dq_i = dqR_2(dq_iph, dq_ip3h);
		double rqR_iph = q_ip1 - dq_i*0.5;
		return rqR_iph;
	}

	//calculating Roe fluxes at each cell boundary (not including ghost cells)
	private double[] calc_flux(double[] q, double[] x, int Nx, int Ng) {
		double[] FF = new double[Nx];

		for (int i = Ng - 1; i < Nx - Ng; i++) {
			double rqL_iph = recos_qL(q[i-1], q[i], q[i+1], x[i-1], x[i], x[i+1]);
			double rqR_iph = recos_qR(q[i], q[i+1], q[i+2], x[i], x[i+1], x[i+2]);

			double rq_iph = 0.5*(rqL_iph + rqR_iph);

			double A = calc_A(rq_iph);
			double eta = calc_eta(A);
			double lambda = calc_lambda(A);
			double omega = calc_omega(rqL_iph, rqR_iph, eta);

			double fR_iph = calc_f(rqR_iph);
			double fL_iph = calc_f(rqL_iph);
			double FF_iph = calc_FF(fR_iph, fL_iph, eta, lambda, omega);
			FF[i] = FF_iph;
		}

		return FF;
	}

	//updating
	private double[] update_q(double[] qn, double[] x, double[] FF, double delta_t, int Nx, int Ng) {
		double dx = x[1] - x[0];
		double[] qnp1 = new double[Nx];

		for (int i = Ng; i < Nx - Ng; i++) {
			qnp1[i] = qn[i] - delta_t/dx*(FF[i] - FF[i-1]);
		}
		
		//constant extrapolation
		for (int i = 0; i < Ng; i++) {
			qnp1[i] = qnp1[Ng];
			qnp1[Nx - i - 1] = qnp1[Nx - Ng - 1];
		}
		
		/*Either zero or first order extrapolation must be commented out
		//first order extrapolation
		for (int i = 0; i < Ng; i++) {
			qnp1[i] = qnp1[Ng]  - (Ng - i)*(qnp1[Ng + 1] - qnp1[Ng]);
			qnp1[Nx - i - 1] = qnp1[Nx - Ng - 1] + (Ng - i)*(qnp1[Nx - Ng - 1] - qnp1[Nx - Ng - 2]);
		}
		*/
		
		return qnp1;
	}

	//moving forward one time step
	private double[] step(double[] q) {

		//half advanced time step
		double[] FF = calc_flux(q, x, Nx, Ng);
		double[] q_nph = update_q(q, x, FF, delta_t/2, Nx, Ng);

		//full advanced time step
		FF = calc_flux(q_nph, x, Nx, Ng);
		double[] q_np1 = update_q(q_nph, x, FF, delta_t, Nx, Ng);

		return q_np1;

	}

	//animates evolution until time tmax
	public void animate(double tmax) {
		StdDraw.setXscale(xmin, xmax);
		StdDraw.setYscale(0.0, 1.5);
		init_fluid();
		numsteps = (int) Math.floor(tmax/delta_t);
		double[] curr = q;
		for (int i = 0; i < numsteps; i++) {
			StdDraw.clear();
			for (int j = 0; j < Nx - 1; j++) {
				StdDraw.line(x[j], curr[j], x[j+1], curr[j+1]);
			}

			//drawing axes
			StdDraw.line(xmin, 0, xmax, 0);
			StdDraw.line(0, 0, 0, 1);
			for (int k = (int) xmin; k <= (int) xmax; k++ ) {
				if (k == 0) {}
				else {
					StdDraw.line((double) k, 0, (double) k, 0.01);
					String I = "" + k;
					StdDraw.text((double) k, 0.025, I);
				}

			}
			for (int k = 1; k <= 5; k++) {

				StdDraw.line(0, 0.2*k, 0.1, 0.2*k);
				String I = String.format("%.1g%n", 0.2*k);
				StdDraw.text(0.35, 0.2*k, I);
			}
			StdDraw.show(30);
			curr = step(curr);
		}
	}

	//displays the solution at a time tmax
	public void time(double tmax) {
		init_fluid();
		numsteps = (int) Math.floor(tmax/delta_t);
		double[] curr = q;
		for (int i = 0; i < numsteps; i++) {
			curr = step(curr);
		}
		StdDraw.setXscale(xmin, xmax);
		for (int j = 0; j < Nx - 1; j++) {
			StdDraw.line(x[j], curr[j], x[j+1], curr[j+1]);
		}
		//drawing axes
		StdDraw.line(xmin, 0, xmax, 0);
		StdDraw.line(0, 0, 0, 1);
		for (int i = (int) xmin; i <= (int) xmax; i++ ) {
			if (i == 0) {}
			else {
				StdDraw.line((double) i, 0, (double) i, 0.01);
				String I = "" + i;
				StdDraw.text((double) i, 0.025, I);
			}

		}
		for (int i = 1; i <= 5; i++) {

			StdDraw.line(0, 0.2*i, 0.1, 0.2*i);
			String I = String.format("%.1g%n", 0.2*i);
			StdDraw.text(0.35, 0.2*i, I);
		}
	}

	//returns values of q at each mesh point after time tmax
	public double[] returnq(double tmax) {
		init_fluid();
		numsteps = (int) Math.floor(tmax/delta_t);
		double[] curr = q;
		for (int i = 0; i < numsteps; i++) {
			curr = step(curr);
		}
		return curr;
	}

	//graphs the convergence factor as a function of x
	public static void convergence(int Nx, double time) {
		double dt = 5.0/(Nx * 20);
		Burgers one = new Burgers(Nx, dt, true);
		Burgers two = new Burgers(2*Nx, dt, true);
		Burgers four = new Burgers(4*Nx, dt, true);

		double[] One = one.returnq(time);
		double[] Two = two.returnq(time);
		double[] Four = four.returnq(time);
		double[] Q = new double[Nx];
		StdDraw.setXscale(one.xmin, one.xmax);
		StdDraw.setYscale(0, 10);
		//drawing axes
		StdDraw.line(one.xmin, 0, one.xmax, 0);
		StdDraw.line(0, 0, 0, 10);
		for (int i = (int) one.xmin; i <= (int) one.xmax; i++ ) {
			if (i == 0) {}
			else {
				StdDraw.line((double) i, 0, (double) i, 0.1);
				String I = "" + i;
				StdDraw.text((double) i, 0.25, I);
			}

		}
		for (int i = 0; i <= 10; i++) {
			if (i == 0) {}
			else {
				StdDraw.line(0, (double) i, 0.1, (double) i);
				String I = "" + i;
				StdDraw.text(0.25, (double) i, I);
			}
		}

		StdDraw.setPenRadius(0.005);
		for (int i = 0; i < Nx; i++) {
			double e1 = Math.abs(One[i] - ((Two[2*i] + Two[(2*i)+1])/2));
			double e2 = Math.abs((((Two[2*i] + Two[(2*i)+1])/2) - (Four[4*i] + Four[4*i + 1] + Four[4*i + 2] + Four[4*i + 3])/4));
			Q[i] = e1/e2;
			StdDraw.point(one.x[i], e1/e2);
		}

	}

	public static void main(String args[]){

		//Burgers(double xmin, double xmax, int Nx, double dt, boolean smooth)
		//Burgers hi = new Burgers(-2, 2, 200, .01, true);
		Burgers hi = new Burgers(-2, 2, 400, .01, false);
		hi.animate(2);

	}
}
