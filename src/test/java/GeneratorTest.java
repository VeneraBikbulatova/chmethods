import org.example.Generator;
import org.junit.Test;

public class GeneratorTest {
    @Test
    public void Test1() {
        int N = 3;
        double ALPHA = 1.0;
        double BETA = 1000.0;
        int n = 3;
        double[][] a = new double[N][];

        for(int i = 0; i < n; ++i) {
            a[i] = new double[n+1];
        }

        double[][] a_inv = new double[n][];

        for(int i = 0; i < n; ++i) {
            a_inv[i] = new double[n];
        }

        Generator g = new Generator();
        double[] values = g.myGen(a, a_inv, n, ALPHA, BETA, 1, 2, 1, 1);
        g.printMatrix(a, n);
        g.printMatrix(a_inv, n);
        g.solve(a, n);
        g.printTable(values);
    }

    @Test
    public void test2(){
        double[][] a = {{1, -2, 1, 6}, {2, -1, 3, 11}, {2, 3, -4, -7}};
        Generator g = new Generator();
        g.solve(a,3);
    }
}

