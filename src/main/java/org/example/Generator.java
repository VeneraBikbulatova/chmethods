package org.example;

public class Generator {
    public Generator() {
    }

    public void printMatrix(double[][] a, int n) {
        System.out.println(" ");
        //  System.out.println(" a:");
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

    public double vectorInfNorm(double[] expect, int n) {
        double max = 0;
        for (int k = 0; k < n; k++) {
            if (max < Math.abs(expect[k])) {
                max = Math.abs(expect[k]);
            }
        }
        return max;
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

    public double[] myGen(double[][] a, double[][] a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant, int schema) {
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


        double[] expect = new double[n];
        for (int k = 0; k < n; k++) {
            expect[k] = 1;
        }
//        double[] result = solve(returnA_B(a, expect, n), n);
        double[] result = gaussZeidel(returnA_B(a, expect, n), n);
        double max = -1;
        double max_e = -1;
        for (int k = 0; k < n; k++) {
            if (max < Math.abs(expect[k] - result[k])) {
                max = Math.abs(expect[k] - result[k]);
            }
            if (max_e < Math.abs(expect[k])) {
                max_e = Math.abs(expect[k]);
            }
        }

        double[] values = new double[9];
        values[0] = alpha;
        values[1] = beta;
        values[2] = matrixInfNorm(a, n);
        values[3] = matrixInfNorm(a_inv, n);
        values[4] = matrixInfNorm(a, n) * matrixInfNorm(a_inv, n);
        values[5] = max;
        values[6] = (max / max_e);
        double[] r_table = new double[n];
        double mul = 0;
        double[] b = new double[n];
        double sum;
        for (int k = 0; k < n; k++) {
            sum = 0;
            for (int l = 0; l < n; l++) {
                sum += a[k][l] * expect[l];
            }
            b[k] = sum;
        }
        for (int k = 0; k < n; k++) {
            mul = 0;
            for (int l = 0; l < n; l++) {
                mul += a[k][l] * result[l];
            }
            r_table[k] = mul - b[k];
        }
        values[7] = vectorInfNorm(r_table, n);
        values[8] = vectorInfNorm(r_table, n) / vectorInfNorm(b, n);
        return values;
    }

    public double[][] returnA_B(double[][] a, double[] res, int n) {
        double[][] a_b = new double[n][n + 1];
        double sum;
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                sum += a[i][j] * res[j];
                a_b[i][j] = a[i][j];
            }
            a_b[i][n] = sum;
        }
        return a_b;
    }

    //2 лаба - метод гаусса-зейделя
    public double[] gaussZeidel(double[][] a_b, int n) {
        double[] solved = new double[n];
        double[] current = new double[n];
        double[] difference = new double[n];
        double min;
        for (int i = 0; i < n; i++) {
            solved[i] = 1.;
        }

        do {
            min = 100;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    current[i] += a_b[i][j] * solved[j];
                }
                current[i] = solved[i] / a_b[i][0];
                difference[i] = Math.abs(current[i]-solved[i]);
                solved[i] = current[i];
                if(solved[i] < min){
                    min = solved[i];
                }
            }
        } while (min > 1.);

        return solved;
    }

    //Метод отражений
    public double[] solve(double[][] a_b, int n) {
        double[] solved = new double[n];
        double[] s = new double[n];
        int count;
        double[] w = new double[n];
        double s_i;
        double ro;
        double alpha;
        int[] i = new int[n];
        double[][] u = new double[n][n];
//        //Создали единичную матрицу
        for (count = 0; count < n; count++) {
            //Ввели s
            for (int j = 0; j < n; j++) {
                s[j] = a_b[j][count];
            }

            /* ================================================= */
            for (int j = 0; j < count; j++) {
                s[j] = 0.0;
            }
            /* ================================================= */


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
            for (int j = 0; j < count; j++) {
                s[j] = 0;
            }
            s_i = s[count];

            ro = Math.sqrt(2.0 * alpha * alpha + 2.0 * Math.abs(alpha) * Math.abs(s_i));
            for (int j = 0; j < n; j++) {
                w[j] = (1.0 / ro) * (s[j] - alpha * i[j]);
            }

            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    if (k == j) {
                        u[k][j] = 1.0 - 2.0 * w[k] * w[j];
                    } else {
                        u[k][j] = -2.0 * w[k] * w[j];
                    }
                }
            }

            a_b = matrixMul34(u, a_b, n, count);

        }
        double sum;
        for (int j = n - 1; j >= 0; j--) {
            sum = 0;
            for (int k = j + 1; k < n; k++) {
                sum += a_b[j][k] * solved[k];
            }
            solved[j] = (a_b[j][n] - sum) / a_b[j][j];
        }
        return solved;
    }

    public void printMatrix34(double[][] a, int n) {
        System.out.println(" ");

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n + 1; ++j) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }
    }

    public void printTable(double[] values) {
        System.out.println("| alpha ||   beta   ||        nrmA         ||        nrmA_        ||         nu          ||          z          ||         ksi          ||          r           ||         ro         |");
        System.out.println("|  " + values[0] + "  ||  " + values[1] + "  ||  " + values[2] + "  ||  " + values[3] + "  ||  " + values[4] + "  ||  " + values[5]
                + "  ||  " + values[6] + "  ||  " + values[7] + "  ||  " + values[8] + "|");
    }

    public double[][] matrixMul34(double[][] a, double[][] b, int n, int count) {
        double[][] c = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                if (i < count && j < count) {
                    c[i][j] = b[i][j];
                } else {
                    c[i][j] = 0.0;
                    for (int k = 0; k < n; k++) {
                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        }
        return c;
    }

}
