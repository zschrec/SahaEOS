

public class Saha {


	static double tol = 1.0e-4; 

	static double m_e = 9.10938291e-28;
	static double h = 6.62606957e-27;
	static double rho = 1.0e-9;
	static double k = 1.3806488e-16;
	static double k_eV = 8.6173324e-5;
	static double N_o = 6.0221413e23;
	static double pi = Math.PI;

	static double CONST = Math.pow(2.0*pi*m_e*k/Math.pow(h,2.0), 3.0/2.0);

	static double[] x = { 1.0 };
	static double[] A = { 1.0080 };
	static double[] m = { 1.673e-24 };
	static double[] n = { N_o*rho*x[0]/A[0] };

	static double[] X_h_eV = { 13.5984 };
	static double[] X_h = { 2.1787e-11 };
	static double[] B_h = { Math.pow(10.0, 0.3), Math.pow(10.0, 0.0) };

	static double N_e = 1.0e5;
	static double N_e1 = 0.0;

	static double N_e_initial, N_e_last;
	
	static double T_N_e, T_N_e_last; 	
	static double g_N_e, g_N_e_last; 	
	static double G_N_e, G_N_e_last; 	

	public static void main(String[] args) {

		N_e_initial = N_e;
		N_e_last = N_e;

		for(int p = 10; p <= 3500; p += 10) {

			N_e = N_e_initial;
			N_e_last = 100.0;
			N_e1 = 0.0;

			for(;;) {
				T_N_e = T(N_e, 1.0*p);
				T_N_e_last = T(N_e_last, 1.0*p);

				g_N_e = g(N_e, T_N_e);
				g_N_e_last = g(N_e_last, T_N_e_last);
				
				G_N_e = G(N_e, g_N_e);
				G_N_e_last = G(N_e_last, g_N_e_last);


				System.out.println("P " + p);
				System.out.println("N_e " + N_e);
				System.out.println("N_e_last " + N_e_last);
				System.out.println("T_N_e " + T_N_e);
				System.out.println("T_N_e_last " + T_N_e_last);
				System.out.println("g_N_e " + g_N_e);
				System.out.println("g_N_e_last " + g_N_e_last);
				System.out.println("G_N_e " + G_N_e);
				System.out.println("G_N_e_last " + G_N_e_last);

				N_e1 = N_e - ( (G_N_e * ( N_e - N_e_last ) ) / (G_N_e - G_N_e_last) );

				if((Math.abs(N_e - N_e1)/N_e) < tol) {
					break;
				}
	
				if( Math.abs(N_e_last - N_e1) / N_e_last > 1.0e-3 ) {

					//System.out.println("HERE " + Math.abs(N_e_last - N_e1) / N_e_last + " " +
					//	Math.abs(N_e - N_e1) / N_e);
					N_e_last = N_e;
					N_e = N_e1;

				} else {

					//System.out.println("THERE " + Math.abs(N_e_last - N_e1) / N_e_last + " " +
					//	Math.abs(N_e - N_e1) / N_e);
					
					N_e_last = N_e;
					N_e = (N_e + N_e1) / 2.0;
				
				}

			
			}
			
			N_e = N_e1;
			T_N_e = T(N_e, 1.0*p);
			g_N_e = g(N_e, T_N_e);
			G_N_e = 1.0 / ((1.0/g_N_e) + 1.0) ;

			//System.out.println(T_N_e + " " + G_N_e + " " + p);
			

		}


	}

	static double g(double N_e, double T_N_e) {
	
		return CONST * (2.0 * B_h[1] / B_h[0]) * Math.pow(T_N_e, 3.0/2.0) * 
			Math.exp(-X_h_eV[0]/(k_eV*T_N_e)) * (1.0/N_e);
	
	}

	static double G(double N_e, double g) {

		return N_e - (n[0] / ((1.0/g) + 1.0));

	}

	static double T(double N_e, double P) {

		return P / (k * (n[0] + N_e));

	}

}
