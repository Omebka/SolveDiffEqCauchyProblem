import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;

public class SolveDiffEqCauchyProblem {
    private double xMax;
    private double h;
    private int precision;

    /**
     * Объект для решения дифференциального уравнения.
     *
     * @param xMax правая граница диапазона значений аргумента.
     * @param h шаг сетки значений аргумента.
     * @param precision точность вычислений (количество знаков после запятой у значений).
     */
    public SolveDiffEqCauchyProblem(double xMax, double h, int precision) {
        this.xMax = xMax;
        this.h = h;
        this.precision = precision;
    }

    public double getXMax() {
        return xMax;
    }

    public void setXMax(double xMax) {
        this.xMax = xMax;
    }

    public double getH() {
        return h;
    }

    public void setH(double h) {
        this.h = h;
    }

    public int getPrecision() {
        return precision;
    }

    public void setPrecision(int precision) {
        this.precision = precision;
    }

    /**
     * Метод вычисления значений узлов сетки аргумента.
     *
     * @param linDiffEqCauchyProblem задача Коши с линейным дифференциальным уравнением вида
     *                               a_0 * y^(n) + a_1 * y^(n - 1) + a_2 * y^(n - 2) + ... +
     *                               a_(n - 2) * y'' + a_(n - 1) * y' + a_n * y = 0
     *                               (без слагаемого с 'x' и 'f(x)').
     *
     * @return Значения узлов сетки аргумента.
     */
    public double[] calcX(LinDiffEqCauchyProblem linDiffEqCauchyProblem) {
        double xMin = linDiffEqCauchyProblem.getXMin();
        int numberOfSteps = calcNumberOfSteps(xMin, xMax, h);

        double[] x = new double[numberOfSteps];
        double x0 = round(xMin);

        for (int i = 0; i < numberOfSteps; i++) {
            x[i] = round(x0);
            x0 += round(h);
        }

        return x;
    }

    /**
     * Метод решения задачи Коши с линейным дифференциальным уравнением n-го порядка вида
     * a_0 * y^(n) + a_1 * y^(n - 1) + a_2 * y^(n - 2) + ... + a_(n - 2) * y'' + a_(n - 1) * y' + a_n * y = 0
     * (без слагаемого с 'x' и 'f(x)')
     * методом Рунге-Кутта 4-го порядка.
     *
     * @param linDiffEqCauchyProblem задача Коши с линейным дифференциальным уравнением вида
     *                               a_0 * y^(n) + a_1 * y^(n - 1) + a_2 * y^(n - 2) + ... +
     *                               a_(n - 2) * y'' + a_(n - 1) * y' + a_n * y = 0
     *                               (без слагаемого с 'x' и 'f(x)').
     *
     * @return Искомые значения функции 'y(x)' в узлах сетки.
     */
    public double[] solveLinDiffEqRungeKutta(LinDiffEqCauchyProblem linDiffEqCauchyProblem) {
        double[] coefficients = linDiffEqCauchyProblem.getCoefficients();
        double[] initialValues = linDiffEqCauchyProblem.getInitialValues();
        double xMin = linDiffEqCauchyProblem.getXMin();

        int stepsRungeKutta = 4;
        double[] y = null;

        if (coefficients.length != initialValues.length + 1 ||
                xMin > xMax ||
                h <= 0) {
            System.out.println("Unsolvable, because of entering incorrect data.");
        } else {
            int numberOfSteps = calcNumberOfSteps(xMin, xMax, h);
            y = new double[numberOfSteps];

            y[0] = round(initialValues[0]);

            // Проверка максимального порядка производных.
            // От этого зависит, нужно ли вводить новые функции: если макс. порядок - 1, то не нужно. Иначе, нужно.
            if (coefficients.length < 3) {
                for (int i = 1; i < numberOfSteps; i++) {
                    // k - массив значений, постепенно приближающих значения главной функции 'y' к искомым,
                    // полученных на 4-х этапах метода Рунге-Кутта.
                    double[] k = new double[stepsRungeKutta];

                    // 1-й шаг метода.
                    k[0] = round(h * fLinearWithoutX(y[i - 1], coefficients));

                    // 2-й и 3-й шаги метода.
                    for (int j = 0; j < 2; j++) {
                        double yIMinus1AddKDividedTo2 = y[i - 1] + k[j] / 2.0;

                        k[j + 1] = round(h * fLinearWithoutX(yIMinus1AddKDividedTo2, coefficients));
                    }

                    // 4-й шаг метода.
                    double yIMinus1AddK3 = y[i - 1] + k[2];

                    k[3] = round(h * fLinearWithoutX(yIMinus1AddK3, coefficients));

                    // Получаем искомое значение главной функции.
                    y[i] = round(y[i - 1] + deltaRungeKutta(k));
                }
            } else {
                // Вводим новые функции z_i.
                // z - Массив значений введенных функций.
                double[] z = new double[initialValues.length - 1];

                for (int i = 0; i < z.length; i++) {
                    z[i] = round(initialValues[i + 1]);
                }

                for (int i = 1; i < numberOfSteps; i++) {
                    // k - массив массивов значений,
                    // постепенно приближающих значения главной функции 'y' и введенных функций 'z_i' к искомым,
                    // полученных на 4-х этапах метода Рунге-Кутта.
                    double[][] k = new double[coefficients.length - 1][stepsRungeKutta];


                    // 1-й шаг метода.
                    k[0][0] = round(h * fLinearWithoutX(y[i - 1], z, coefficients));

                    for (int j = 1; j < k.length; j++) {
                        k[j][0] = round(h * z[j - 1]);
                    }


                    // 2-й и 3-й шаги метода.
                    for (int j = 0; j < 2; j++) {
                        double yIMinus1AddKDividedTo2 = y[i - 1] + k[1][j] / 2.0;

                        // z1 - вспомогательный массив значений введенных функций z_i, полученных на 2-м и 3-м этапах.
                        double[] z1 = new double[z.length];

                        for (int q = 0; q < z1.length; q++) {
                            if (q == z1.length - 1) {
                                z1[q] = z[q] + k[0][j] / 2.0;
                            } else {
                                z1[q] = z[q] + k[q + 2][j] / 2.0;
                            }
                        }

                        k[0][j + 1] = round(h * fLinearWithoutX(yIMinus1AddKDividedTo2, z1, coefficients));

                        for (int q = 1; q < k.length; q++) {
                            k[q][j + 1] = round(h * z[q - 1]);
                        }
                    }


                    // 4-й шаг метода.
                    double yIMinus1AddK3 = y[i - 1] + k[1][2];

                    // z3 - вспомогательный массив значений введенных функций z_i, полученных на 4-м этапе.
                    double[] z3 = new double[z.length];

                    for (int j = 0; j < z3.length; j++) {
                        if (j == z3.length - 1) {
                            z3[j] = z[j] + k[0][2];
                        } else {
                            z3[j] = z[j] + k[j + 2][2];
                        }
                    }

                    k[0][3] = round(h * fLinearWithoutX(yIMinus1AddK3, z3, coefficients));

                    for (int j = 1; j < k.length; j++) {
                        k[j][3] = round(h * z[j - 1]);
                    }


                    // Получаем искомое значение главной функции.
                    y[i] = round(y[i - 1] + deltaRungeKutta(k[1]));

                    // Получаем искомые значения введенных функций.
                    for (int j = 0; j < z.length; j++) {
                        if (j == z.length - 1) {
                            z[j] = round(z[j] + deltaRungeKutta(k[0]));
                        } else {
                            z[j] = round(z[j] + deltaRungeKutta(k[j + 2]));
                        }
                    }
                }
            }
        }

        return y;
    }

    /**
     * Метод подсчета количества узлов сетки значений аргумента.
     *
     * @param xMin левая граница диапазона значений аргумента.
     * @param xMax правая граница диапазона значений аргумента.
     * @param h шаг сетки значений аргумента.
     *
     * @return Значение количества узлов сетки значений аргумента.
     */
    private static int calcNumberOfSteps(double xMin, double xMax, double h) {
        return (int) ((xMax - xMin) / h) + 1;
    }

    /**
     * Метод вычисления значения функции,
     * полученной выделением слагаемого с максимальным порядком производной из линейного дифференциального уравнения:
     * переносом всего остального в другую часть уравнения и
     * делением этого всего на коэффициент перед производной максимального порядка.
     *
     * @param y предыдущее значение искомой функции.
     * @param z массив предыдущих значений введенных функций.
     * @param coefficients значения коэффициентов 'a_i' перед соответствующими производными искомой функции 'y'
     *                     в порядке убывания порядка производной (начиная с наивысшего порядка 'n'):
     *                     a_0 * y^(n),
     *                     a_1 * y^(n - 1),
     *                     a_2 * y^(n - 2),
     *                     ...,
     *                     a_(n - 2) * y'',
     *                     a_(n - 1) * y',
     *                     a_n * y.
     *
     * @return Подсчитанное значение искомой функции.
     */
    private static double fLinearWithoutX(double y, double[] z, double[] coefficients) {
        double mainCoef = coefficients[0];
        double result = 0;

        for (int i = 0; i < z.length + 1; i++) {
            if (i == z.length) {
                result += -coefficients[i + 1] / mainCoef * y;
            } else {
                result += -coefficients[i + 1] / mainCoef * z[z.length - 1 - i];
            }
        }

        return result;
    }

    /**
     * Метод вычисления значения функции, полученной выделением слагаемого с производной 1-го порядка
     * из линейного дифференциального уравнения 1-го порядка:
     * переносом всего остального в другую часть уравнения и
     * делением этого всего на коэффициент перед производной 1-го порядка.
     *
     * @param y предыдущее значение искомой функции.
     * @param coefficients значения коэффициентов 'a_i' перед соответствующими производными искомой функции 'y'
     *                     в порядке убывания порядка производной (начиная с наивысшего порядка 1):
     *                     a_0 * y',
     *                     a_1 * y.
     *
     * @return Подсчитанное значение искомой функции.
     */
    private static double fLinearWithoutX(double y, double[] coefficients) {
        return -coefficients[1] / coefficients[0] * y;
    }

    /**
     * Метод подсчета значения приращения между текущим и следующим значениями функции методом Рунге-Кутта.
     *
     * @param k массив значений, постепенно приближающих значения функции к искомым,
     *          полученных на 4-х этапах метода Рунге-Кутта.
     *
     * @return Значение приращения между текущим и следующим значениями функции.
     */
    private static double deltaRungeKutta(double[] k) {
        return (k[0] + 2 * k[1] + 2 * k[2] + k[3]) / 6.0;
    }

    /**
     * Метод округления значения 'value' до определенного количества знаков после запятой.
     *
     * @param value число, которое нужно округлить.
     *
     * @return Округленное значение переменной 'value'.
     */
    private double round(double value) {
        double powOf10 = Math.pow(10, precision);
        return (double) Math.round(value * powOf10) / powOf10;
    }

    /**
     * Метод построения графика функции.
     *
     * @param x массив double значений аргумента.
     * @param y массив double значений функции.
     * @param chartName название окна графика.
     */
    public static void graphics(double[] x, double[] y, String chartName) {
        XYSeries seriesY = new XYSeries("y(x)");

        for (int i = 0; i < y.length; i++) {
            seriesY.add(x[i], y[i]);
        }

        XYSeriesCollection xy = new XYSeriesCollection(seriesY);
        createChart(xy, chartName);
    }

    /**
     * Метод создания окна графика функции.
     *
     * @param xy коллекция значений аргументов 'x' и их значений функции 'y'.
     * @param chartName название окна графика.
     */
    private static void createChart(XYSeriesCollection xy, String chartName) {
        JFreeChart chart = ChartFactory
                .createXYLineChart("", "x", "", xy, PlotOrientation.VERTICAL,
                        true, true, true);
        JFrame frame = new JFrame(chartName);
        frame.getContentPane().add(new ChartPanel(chart));
        frame.setSize(400, 300);
        frame.setVisible(true);

        frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
    }
}
