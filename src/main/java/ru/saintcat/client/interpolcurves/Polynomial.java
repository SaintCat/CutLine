/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ru.saintcat.client.interpolcurves;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 *
 * @author Chernyshov
 */
public class Polynomial {

    public DenseMatrix64F Ux;
    public DenseMatrix64F Uy;
    public DenseMatrix64F Sx;
    public DenseMatrix64F Sy;

    public List<Vector2D> Normal(List<Vector2D> points, double dt, boolean secondDerivative, boolean loop) {
        List<Vector2D> normalPoints = new ArrayList<>(points);
        if (loop) {
            normalPoints.add(normalPoints.get(0));
        }
        int n = normalPoints.size();
        normalPoints.add(new Vector2D(0.0, 0.0));
        if (secondDerivative) {
            normalPoints.add(new Vector2D(0.0, 0.0));
        }

        int m = normalPoints.size() - 1;
        double[] tmpX = new double[normalPoints.size()];
        double[] tmpY = new double[normalPoints.size()];
        for (int i = 0; i < normalPoints.size(); i++) {
            tmpX[i] = normalPoints.get(i).x;
            tmpY[i] = normalPoints.get(i).y;
        }

        Ux = DenseMatrix64F.wrap(1, normalPoints.size(), tmpX);
        Uy = DenseMatrix64F.wrap(1, normalPoints.size(), tmpY);

        System.out.println("Ux");
        Ux.print();
        System.out.println("Uy");
        Uy.print();
        DenseMatrix64F Qmn = new DenseMatrix64F(m + 1, m + 1);

        if (secondDerivative) {
            Qmn = GetQmndd(m, n);
        } else {
            Qmn = GetQmnd(m, n);
        }

        System.out.println("Qmn");
        Qmn.print();
        if (CommonOps.det(Qmn) == 0) {
            return new ArrayList<Vector2D>();
        }

        CommonOps.invert(Qmn);
        Qmn.print();
        Sx = new DenseMatrix64F(Ux.numRows, Qmn.numCols);
        Sy = new DenseMatrix64F(Ux.numRows, Qmn.numCols);
        CommonOps.mult(Ux, Qmn, Sx);
        CommonOps.mult(Uy, Qmn, Sy);

        List<Vector2D> res = new ArrayList<Vector2D>();
        for (double t = 0.0; t <= n - 1; t += dt) {
            double sumX = 0.0;
            double sumY = 0.0;
            double[] Tm = GetTm(m, t);

            for (int i = 0; i <= m; i++) {
                sumX += Sx.get(0, i) * Tm[i];
                sumY += Sy.get(0, i) * Tm[i];
            }

            res.add(new Vector2D(sumX, sumY));
        }
        return res;
    }

    public List<Vector2D> Chord(List<Vector2D> points, double dt, boolean secondDerivative, boolean loop) {
        List<Vector2D> normalPoints = new ArrayList<Vector2D>(points);
        if (loop) {
            normalPoints.add(normalPoints.get(0));
        }
        int n = normalPoints.size();
        normalPoints.add(new Vector2D(0.0, 0.0));
        if (secondDerivative) {
            normalPoints.add(new Vector2D(0.0, 0.0));
        }

        int m = normalPoints.size() - 1;
        double[] tmpX = new double[normalPoints.size()];
        double[] tmpY = new double[normalPoints.size()];
        for (int i = 0; i < normalPoints.size(); i++) {
            tmpX[i] = normalPoints.get(i).x;
            tmpY[i] = normalPoints.get(i).y;
        }

        Ux = DenseMatrix64F.wrap(1, normalPoints.size(), tmpX);
        Uy = DenseMatrix64F.wrap(1, normalPoints.size(), tmpY);

        DenseMatrix64F Qmn = new DenseMatrix64F(m + 1, m + 1);

        if (secondDerivative) {
            Qmn = GetChordQmndd(m, normalPoints, n);
        } else {
            Qmn = GetChordQmnd(m, normalPoints, n);
        }
        if (CommonOps.det(Qmn) == 0) {
            return new ArrayList<Vector2D>();
        }

        CommonOps.invert(Qmn);

        Sx = new DenseMatrix64F(Ux.numRows, Qmn.numCols);
        Sy = new DenseMatrix64F(Ux.numRows, Qmn.numCols);
        CommonOps.mult(Ux, Qmn, Sx);
        CommonOps.mult(Uy, Qmn, Sy);

        List<Vector2D> res = new ArrayList<Vector2D>();

        double length = 0;
        for (int i = 0; i < n - 1; i++) {
            length += Math.abs(getLength(normalPoints.get(i + 1).x, normalPoints.get(i + 1).y, normalPoints.get(i).x, normalPoints.get(i).y));
        }
        for (double t = 0; t <= length; t += length / 10000) {
            double sumX = 0.0;
            double sumY = 0.0;
            double[] Tm = GetTm(m, t);

            for (int j = 0; j <= m; j++) {
                sumX += Sx.get(0, j) * Tm[j];
                sumY += Sy.get(0, j) * Tm[j];
            }

            res.add(new Vector2D(sumX, sumY));
        }
        return res;
    }

    private static double getLength(double a, double a1, double b, double b1) {
        return Math.sqrt(Math.pow(b1 - a1, 2) + Math.pow(b - a, 2));
    }

    public DenseMatrix64F GetQmn(int m) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        for (int i = 0; i <= m; i++) {
            double[] tm = GetTm(m, i);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }

        }

        return T;
    }

    public DenseMatrix64F GetQmnd(int m, int n) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        for (int i = 0; i < m; i++) {
            double[] tm = GetTm(m, i);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }
        }
        System.out.println("T");
        T.print();
        // first derivative
        double[] resT = new double[m + 1];
        double[] tmpT0 = GetTmd(m, 0);
        System.out.println("tmpT0");
        System.out.println(Arrays.toString(tmpT0));
        System.out.println("tmpTn");
        double[] tmpTn = GetTmd(m, n - 1);
        System.out.println(Arrays.toString(tmpTn));
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m, resT[j]);
        }
        return T;
    }

    public DenseMatrix64F GetQmndd(int m, int n) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        for (int i = 0; i < m - 1; i++) {
            double[] tm = GetTm(m, i);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }
        }
        // first derivative
        double[] resT = new double[m + 1];
        double[] tmpT0 = GetTmd(m, 0);
        double[] tmpTn = GetTmd(m, n - 1);
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m - 1, resT[j]);
        }
        // second derivative
        resT = new double[m + 1];
        tmpT0 = GetTmdd(m, 0);
        tmpTn = GetTmdd(m, n - 1);
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m, resT[j]);
        }
        return T;
    }

    public double[] GetTm(int m, double t) {
        double[] tmp = new double[m + 1];
        for (int i = 0; i <= m; i++) {
            tmp[i] = Math.pow(t, i);
        }
        return tmp;
    }

    public double[] GetTmd(int m, double t) {
        double[] tmp = new double[m + 1];
        for (int i = 0; i <= m; i++) {
            double pow = Math.pow(t, i - 1);
            if (Double.isInfinite(pow)) {
                pow = 0.0;
            }
            tmp[i] = i * pow;
        }
        return tmp;
    }

    public double[] GetTmdd(int m, double t) {
        double[] tmp = new double[m + 1];
        for (int i = 0; i <= m; i++) {
            double pow = Math.pow(t, i - 2);
            if (Double.isInfinite(pow)) {
                pow = 0.0;
            }
            tmp[i] = i * (i - 1) * pow;
        }
        return tmp;
    }

    public DenseMatrix64F GetChordQmn(int m, List<Vector2D> points) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        double length = 0;
        for (int i = 0; i <= m; i++) {
            double[] tm = GetTm(m, length);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }
            length += Math.abs(getLength(points.get(i + 1).x, points.get(i + 1).y, points.get(i).x, points.get(i).y));
        }

        return T;
    }

    public DenseMatrix64F GetChordQmnd(int m, List<Vector2D> points, int n) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        double length = 0;
        double nLength = 0;
        for (int i = 0; i < m; i++) {
            double[] tm = GetTm(m, length);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }
            if (i == n - 1) {
                nLength = length;
            }
            length += Math.abs(getLength(points.get(i + 1).x, points.get(i + 1).y, points.get(i).x, points.get(i).y));
        }
        // first derivative
        double[] resT = new double[m + 1];
        double[] tmpT0 = GetTmd(m, 0);
        length += Math.abs(getLength(points.get(n - 1).x, points.get(n - 1).y, points.get(0).x, points.get(0).y));
        double[] tmpTn = GetTmd(m, nLength);
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m, resT[j]);
        }
        return T;
    }

    public DenseMatrix64F GetChordQmndd(int m, List<Vector2D> points, int n) // Матрица Вандермонда
    {
        DenseMatrix64F T = new DenseMatrix64F(m + 1, m + 1);
        double length = 0;
        double nLength = 0;
        for (int i = 0; i < m; i++) {
            if (i == n - 1) {
                nLength = length;
            }
            double[] tm = GetTm(m, length);
            for (int j = 0; j < tm.length; j++) {
                T.set(j, i, tm[j]);
            }
            length += Math.abs(getLength(points.get(i + 1).x, points.get(i + 1).y, points.get(i).x, points.get(i).y));
        }
        // first derivative
        double[] resT = new double[m + 1];
        double[] tmpT0 = GetTmd(m, 0);
        length += Math.abs(getLength(points.get(n - 1).x, points.get(n - 1).y, points.get(0).x, points.get(0).y));
        double[] tmpTn = GetTmd(m, nLength);
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m - 1, resT[j]);
        }
        // second derivative
        resT = new double[m + 1];
        tmpT0 = GetTmdd(m, 0);
        length += Math.abs(getLength(points.get(n - 1).x, points.get(n - 1).y, points.get(0).x, points.get(0).y));
        tmpTn = GetTmdd(m, nLength);
        for (int i = 0; i < tmpT0.length; i++) {
            resT[i] = tmpTn[i] - tmpT0[i];
        }
        for (int j = 0; j < resT.length; j++) {
            T.set(j, m, resT[j]);
        }
        return T;
    }
}
