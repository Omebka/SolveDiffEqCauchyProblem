public class LinDiffEqCauchyProblem {
    private double xMin;
    private double[] coefficients;
    private double[] initialValues;

    /**
     * Объект - задача Коши с линейным дифференциальным уравнением вида
     * a_0 * y^(n) + a_1 * y^(n - 1) + a_2 * y^(n - 2) + ... + a_(n - 2) * y'' + a_(n - 1) * y' + a_n * y = 0
     * (без слагаемого с 'x' и 'f(x)').
     *
     * @param xMin левая граница диапазона значений аргумента.
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
     *                     Количество элементов - n + 1.
     * @param initialValues значения начальных условий задачи Коши в точке 'xMin'
     *                      в порядке возрастания порядка производной (начиная с наименьшего порядка '0'):
     *                      y(xMin),
     *                      y'(xMin),
     *                      y''(xMin),
     *                      ...,
     *                      y^(n - 2)(xMin),
     *                      y^(n - 1)(xMin).
     *
     *                      Количество элементов - n.
     */
    public LinDiffEqCauchyProblem(double xMin, double[] coefficients, double[] initialValues) {
        this.xMin = xMin;
        this.coefficients = coefficients;
        this.initialValues = initialValues;
    }

    public double getXMin() {
        return xMin;
    }

    public void setXMin(double xMin) {
        this.xMin = xMin;
    }

    public double[] getCoefficients() {
        return coefficients;
    }

    public void setCoefficients(double[] coefficients) {
        this.coefficients = coefficients;
    }

    public double[] getInitialValues() {
        return initialValues;
    }

    public void setInitialValues(double[] initialValues) {
        this.initialValues = initialValues;
    }
}
