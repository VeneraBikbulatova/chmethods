package org.example;

public class Generator {
    public Generator() {
    }

    public void printMatrix(double[][] a, int n) {
        System.out.println(" ");
        System.out.println(" a:");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }
    }

    public double matrixInfNorm(double[][] a, int n) {
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j) {
                s += Math.abs(a[i][j]);
            }
            if (s > norm) {
                norm = s;
            }
        }
        return norm;
    }

    public void matrixMul(double[][] a, double[][] b, double[][] c, int n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                c[i][j] = 0.0;
                for (int k = 0; k < n; ++k) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
    }

    //ортогональная матрица
    public void Q_matrix(double[][] Q, int n, int schema) {
        double next = 1.0;
        int i;
        double q;
        for (int j = 0; j < n - 1; ++j) {
            double curr = next++;
            q = 1.0 / Math.sqrt(curr * next);
            for (i = 0; i <= j; ++i) {
                Q[i][j] = q;
            }
            Q[j + 1][j] = -Math.sqrt(curr / next);
            for (i = j + 2; i < n; ++i) {
                Q[i][j] = 0.0;
            }
        }
        q = 1.0 / Math.sqrt(n);
        for (i = 0; i < n; ++i) {
            Q[i][n - 1] = q;
        }
    }

    public void myGen(double[][] a, double[][] a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema) {
        System.out.println("   M A T R I X  G E N.  ");
        System.out.println("              N = " + n);
        System.out.println(" | lambda_min | = " + alpha);
        System.out.println(" | lambda_max | = " + beta);
        double[] lambda = new double[n];
        System.out.println(" sign_law = " + sign_law);
        double[] sign = new double[n];
        int i;
        for (i = 0; i < n; ++i) {
            sign[i] = 1.0;
        }
        label478:
        switch (sign_law) {
            case -1:
                i = 0;

                while (true) {
                    if (i >= n) {
                        break label478;
                    }

                    sign[i] = -1.0;
                    ++i;
                }
            case 0:
                sign[0] = 1.0;
                for (i = 1; i < n; ++i) {
                    sign[i] = -sign[i - 1];
                }
        }
        System.out.println(" lambda_law = " + lambda_law);
        double[] kappa = new double[n];
        for (i = 0; i < n; ++i) {
            kappa[i] = (double) i / (double) (n - 1);
        }
        label458:
        switch (lambda_law) {
            case 1:
                System.out.println(" kappa = sqrt( ) ");
                i = 0;
                while (true) {
                    if (i >= n) {
                        break label458;
                    }

                    kappa[i] = Math.sqrt(kappa[i]);
                    ++i;
                }
            case 2:
                System.out.println(" kappa = sin( ) ");
                double pi_half = Math.acos(-1.0) * 0.5;

                for (i = 0; i < n; ++i) {
                    kappa[i] = Math.sin(pi_half * kappa[i]);
                }
        }
        double[] J = new double[n];
        for (i = 0; i < n; ++i) {
            J[i] = sign[i] * ((1.0 - kappa[i]) * alpha + kappa[i] * beta);
        }
        double[] J_inv = new double[n];
        for (i = 0; i < n; ++i) {
            J_inv[i] = 1.0 / J[i];
        }
        double[][] Q = new double[n][];
        for (i = 0; i < n; ++i) {
            Q[i] = new double[n];
        }
        double s;
        double[] aa = new double[3];
        System.out.println(" variant = " + variant);
        int j;
        label426:
        switch (variant) {
            case 0 -> {
                System.out.println(" simmetric matrix:");
                System.out.println(" schema = " + schema);
                if (schema == 1) {
                    this.Q_matrix(Q, n, schema);
                    a[0][0] = 0.0;
                    int k;
                    for (k = 0; k < n; ++k) {
                        a[0][0] += Q[0][k] * J[k] * Q[0][k];
                    }
                    for (j = 1; j < n; ++j) {
                        a[0][j] = 0.0;
                        for (k = j - 1; k < n; ++k) {
                            a[0][j] += Q[0][k] * J[k] * Q[j][k];
                        }
                        a[j][0] = a[0][j];
                    }
                    for (i = 1; i < n; ++i) {
                        a[i][i] = 0.0;
                        for (k = i - 1; k < n; ++k) {
                            a[i][i] += Q[i][k] * J[k] * Q[i][k];
                        }
                        for (j = i + 1; j < n; ++j) {
                            a[i][j] = 0.0;
                            for (k = j - 1; k < n; ++k) {
                                a[i][j] += Q[i][k] * J[k] * Q[j][k];
                            }
                            a[j][i] = a[i][j];
                        }
                    }
                    a_inv[0][0] = 0.0;
                    for (k = 0; k < n; ++k) {
                        a_inv[0][0] += Q[0][k] * J_inv[k] * Q[0][k];
                    }
                    for (j = 1; j < n; ++j) {
                        a_inv[0][j] = 0.0;
                        for (k = j - 1; k < n; ++k) {
                            a_inv[0][j] += Q[0][k] * J_inv[k] * Q[j][k];
                        }
                        a_inv[j][0] = a_inv[0][j];
                    }
                    for (i = 1; i < n; ++i) {
                        a_inv[i][i] = 0.0;
                        for (k = i - 1; k < n; ++k) {
                            a_inv[i][i] += Q[i][k] * J_inv[k] * Q[i][k];
                        }
                        for (j = i + 1; j < n; ++j) {
                            a_inv[i][j] = 0.0;
                            for (k = j - 1; k < n; ++k) {
                                a_inv[i][j] += Q[i][k] * J_inv[k] * Q[j][k];
                            }
                            a_inv[j][i] = a_inv[i][j];
                        }
                    }
                }
            }
            case 1 -> {
                System.out.println(" simple structure matrix:");
                System.out.println(" schema = " + schema);
                if (schema == 1) {
                    a[0][0] = J[0];
                    a[0][1] = -J[1];
                    for (i = 1; i < n - 1; ++i) {
                        a[i][i - 1] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i][i + 1] = -J[i + 1];
                    }
                    a[n - 1][n - 2] = -J[n - 2];
                    a[n - 1][n - 1] = J[n - 1] + J[n - 1];
                    aa[1] = a[0][0];
                    aa[2] = a[0][1];
                    a[0][0] = aa[1] * (double) n + aa[2] * (double) (n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j < n; ++j) {
                        a[0][j] = s * (double) (n - j);
                    }
                    for (i = 1; i < n - 1; ++i) {
                        aa[0] = a[i][i - 1];
                        aa[1] = a[i][i];
                        aa[2] = a[i][i + 1];
                        for (j = 0; j < i; ++j) {
                            a[i][j] = aa[0] * (double) (n - i + 1) + aa[1] * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        }
                        s = aa[0] + aa[1];
                        a[i][i] = s * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j < n; ++j) {
                            a[i][j] = s * (double) (n - j);
                        }
                    }
                    aa[0] = a[n - 1][n - 2];
                    aa[1] = a[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j < n - 1; ++j) {
                        a[n - 1][j] = s;
                    }
                    a[n - 1][n - 1] = aa[0] + aa[1];
                    a_inv[0][0] = J_inv[0];
                    a_inv[0][1] = -J_inv[1];
                    for (i = 1; i < n - 1; ++i) {
                        a_inv[i][i - 1] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i][i + 1] = -J_inv[i + 1];
                    }
                    a_inv[n - 1][n - 2] = -J_inv[n - 2];
                    a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];
                    aa[1] = a_inv[0][0];
                    aa[2] = a_inv[0][1];
                    a_inv[0][0] = aa[1] * (double) n + aa[2] * (double) (n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j < n; ++j) {
                        a_inv[0][j] = s * (double) (n - j);
                    }
                    for (i = 1; i < n - 1; ++i) {
                        aa[0] = a_inv[i][i - 1];
                        aa[1] = a_inv[i][i];
                        aa[2] = a_inv[i][i + 1];
                        for (j = 0; j < i; ++j) {
                            a_inv[i][j] = aa[0] * (double) (n - i + 1) + aa[1] * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        }
                        s = aa[0] + aa[1];
                        a_inv[i][i] = s * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j < n; ++j) {
                            a_inv[i][j] = s * (double) (n - j);
                        }
                    }
                    aa[0] = a_inv[n - 1][n - 2];
                    aa[1] = a_inv[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j < n - 1; ++j) {
                        a_inv[n - 1][j] = s;
                    }
                    a_inv[n - 1][n - 1] = aa[0] + aa[1];
                }
            }
            case 2 -> {
                System.out.println(" J_2 type matrix: must be n > 2");
                System.out.println(" schema = " + schema);
                if (schema == 1) {
                    a[0][0] = J[0];
                    a[0][1] = 1.0 - J[0];
                    a[1][0] = -J[0];
                    a[1][1] = -1.0 + J[0] + J[0];
                    a[1][2] = -J[2];
                    a[2][1] = -J[0];
                    a[2][2] = J[2] + J[2];
                    if (n > 3) {
                        a[2][3] = -J[3];
                    }
                    for (i = 3; i < n - 1; ++i) {
                        a[i][i - 1] = -J[i - 1];
                        a[i][i] = J[i] + J[i];
                        a[i][i + 1] = -J[i + 1];
                    }
                    if (n > 3) {
                        a[n - 1][n - 2] = -J[n - 2];
                        a[n - 1][n - 1] = J[n - 1] + J[n - 1];
                    }
                    aa[1] = a[0][0];
                    aa[2] = a[0][1];
                    a[0][0] = aa[1] * (double) n + aa[2] * (double) (n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j < n; ++j) {
                        a[0][j] = s * (double) (n - j);
                    }
                    for (i = 1; i < n - 1; ++i) {
                        aa[0] = a[i][i - 1];
                        aa[1] = a[i][i];
                        aa[2] = a[i][i + 1];
                        for (j = 0; j < i; ++j) {
                            a[i][j] = aa[0] * (double) (n - i + 1) + aa[1] * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        }
                        s = aa[0] + aa[1];
                        a[i][i] = s * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j < n; ++j) {
                            a[i][j] = s * (double) (n - j);
                        }
                    }
                    aa[0] = a[n - 1][n - 2];
                    aa[1] = a[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j < n - 1; ++j) {
                        a[n - 1][j] = s;
                    }
                    a[n - 1][n - 1] = aa[0] + aa[1];
                    a_inv[0][0] = J_inv[0];
                    a_inv[0][1] = -J_inv[0] * J_inv[0] - J_inv[0];
                    a_inv[1][0] = -J_inv[0];
                    a_inv[1][1] = J_inv[0] * J_inv[0] + J_inv[0] + J_inv[0];
                    a_inv[1][2] = -J_inv[2];
                    a_inv[2][1] = -J_inv[0];
                    a_inv[2][2] = J_inv[2] + J_inv[2];
                    if (n > 3) {
                        a_inv[2][3] = -J_inv[3];
                    }
                    for (i = 3; i < n - 1; ++i) {
                        a_inv[i][i - 1] = -J_inv[i - 1];
                        a_inv[i][i] = J_inv[i] + J_inv[i];
                        a_inv[i][i + 1] = -J_inv[i + 1];
                    }
                    if (n > 3) {
                        a_inv[n - 1][n - 2] = -J_inv[n - 2];
                        a_inv[n - 1][n - 1] = J_inv[n - 1] + J_inv[n - 1];
                    }
                    aa[1] = a_inv[0][0];
                    aa[2] = a_inv[0][1];
                    a_inv[0][0] = aa[1] * (double) n + aa[2] * (double) (n - 1);
                    s = aa[1] + aa[2];
                    for (j = 1; j < n; ++j) {
                        a_inv[0][j] = s * (double) (n - j);
                    }
                    for (i = 1; i < n - 1; ++i) {
                        aa[0] = a_inv[i][i - 1];
                        aa[1] = a_inv[i][i];
                        aa[2] = a_inv[i][i + 1];
                        for (j = 0; j < i; ++j) {
                            a_inv[i][j] = aa[0] * (double) (n - i + 1) + aa[1] * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        }

                        s = aa[0] + aa[1];
                        a_inv[i][i] = s * (double) (n - i) + aa[2] * (double) (n - i - 1);
                        s += aa[2];
                        for (j = i + 1; j < n; ++j) {
                            a_inv[i][j] = s * (double) (n - j);
                        }
                    }
                    aa[0] = a_inv[n - 1][n - 2];
                    aa[1] = a_inv[n - 1][n - 1];
                    s = aa[0] + aa[0] + aa[1];
                    for (j = 0; j < n - 1; ++j) {
                        a_inv[n - 1][j] = s;
                    }
                    a_inv[n - 1][n - 1] = aa[0] + aa[1];
                }
            }
        }
        s = this.matrixInfNorm(a, n);
        System.out.println(" ||  A  || = " + s);
        double norm_inv = this.matrixInfNorm(a_inv, n);
        System.out.println(" ||A_inv|| = " + norm_inv);
        double obusl = s * norm_inv;
        System.out.println(" obusl = " + obusl);
        double[][] r = new double[n][];
        for (i = 0; i < n; ++i) {
            r[i] = new double[n];
        }
        this.matrixMul(a, a_inv, r, n);
        for (i = 0; i < n; ++i) {
            int var10002 = (int) r[i][i]--;
        }
        s = this.matrixInfNorm(r, n);
        System.out.println(" ||R_gen|| = " + s);
    }

    //Метод отражений
    public void solve(double[][] a, double[] b, int n) {
        double[][] a_b = new double[n][n + 1];
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n + 1; j++) {
                if (j == n) {
                    a_b[k][j] = b[k];
                } else {
                    a_b[k][j] = a[k][j];
                }
            }
        }
        double[] s = new double[n];
        int count = 0;
        double[] w = new double[n];
        double s_i = 0;
        double ro;
        double alpha;
        double[][] w_wt = new double[n][n];
        int[] i = new int[n];
        double[][] u = new double[n][n];
        double[][] a_i = new double[n][n + 1];
//        //Создали единичную матрицу
        for (count = 0; count < n; count++) {
            //Ввели s
            for (int j = 0; j < n; j++) {
                if(count == 0){
                    s[j] = a_b[j][count];
                } else {
                    s[j] = a_i[j][count];
                }
            }

            //Ввели i
            for (int j = 0; j < n; j++) {
                i[j] = 0;
            }
            i[count] = 1;
            alpha = 0;
            for (int j = 0; j < n; j++) {
                alpha += s[j] * s[j];
            }
            alpha = Math.sqrt(alpha) * (-s[count] / Math.abs(s[count]));
            //ro or beta
            for(int j = 0; j < count; j++){
                s[j] = 0;
            }
            s_i = s[count];

            ro = Math.sqrt(2 * alpha * alpha + 2 * Math.abs(alpha) * Math.abs(s_i));
            for (int j = 0; j < n; j++) {
                w[j] = (1 / ro) * (s[j] - alpha * i[j]);
            }
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    w_wt[j][k] = w[j] * w[k];
                }
            }

            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    if (k == j) {
                        u[k][j] = 1 - 2 * w_wt[k][j];
                    } else {
                        u[k][j] = 0 - 2 * w_wt[k][j];
                    }
                }
            }
            printMatrix(u, 3);
            System.out.println("u");
            //до этого считает точно правильно!!!
            double[][] ai = new double[n][n + 1];
//            if (count == 0) {
                matrixMul34(u, a_b, a_i, n, count);
//            } else {
//                matrixMul34(u, a_i, a_i, n, count);
//            }
            System.out.println("ALL MATRICES " + count);
            printMatrix34(a_i, n);
            a_b = a_i;
        }
    }

    public void printMatrix34(double[][] a, int n) {
        System.out.println(" ");
        System.out.println(" a:");
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n + 1; ++j) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }
    }


    public void matrixMul34(double[][] a, double[][] b, double[][] c, int n, int count) {
//        for (int i = count; i < n; i++) {
//            for (int j = count; j < n + 1; j++) {
//                c[i][j] = 0.0;
//            }
//        }
        for (int i = count; i < n; i++) {
            for (int j = count; j < n + 1; j++) {
                c[i][j] = 0.0;
                for (int k = 0; k < n; k++) {
//                    if (i > 0 && j == 0) {
//                        c[i][j] = 0.0;
//                    } else {
                        c[i][j] += a[i][k] * b[k][j];
//                    }
                }
            }
        }
    }

}
