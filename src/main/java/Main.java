public class Main {
    public static void main(String[] args) {
        double xMin = 0;
        double[] coefficients = {1, 15, 90, 270, 405, 243};
        double[] initialValues = {0, 3, -9, -8, 0};

        LinDiffEqCauchyProblem linDiffEqCauchyProblem = new LinDiffEqCauchyProblem(xMin, coefficients, initialValues);

        double xMax = 5;
        double h = 0.01;
        int precision = 4;

        SolveDiffEqCauchyProblem solveDiffEqCauchyProblem = new SolveDiffEqCauchyProblem(xMax, h, precision);

        double[] x = solveDiffEqCauchyProblem.calcX(linDiffEqCauchyProblem);
        double[] y = solveDiffEqCauchyProblem.solveLinDiffEqRungeKutta(linDiffEqCauchyProblem);

        if (y != null) {
            SolveDiffEqCauchyProblem.graphics(x, y, "График");
        }
    }
}
