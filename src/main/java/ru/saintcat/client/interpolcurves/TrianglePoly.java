/*
 * TrianglePoly.java
 * 
 * Created: 24.12.2014
 * 
 * Copyright (c) ExpertPB 2014
 * All information contained herein is, and remains the property of
 * ExpertPB and its suppliers, if any.
 */
package ru.saintcat.client.interpolcurves;

import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Roman Chernyshev
 */
public class TrianglePoly {

    public static List<List<Vector2D>> triIteratePoly(List<Vector2D> P, Vector2D o, Vector2D V) {
        int n = P.size();
        V = V.normal();
        List<List<Vector2D>> result = new ArrayList<>();
        List<Double> L = new ArrayList<>();
        List<Vector2D> Ll = new ArrayList<>();
        List<Vector2D> Lr = new ArrayList<>();
        List<Double> fies = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            fies.add(nfV(o, V, P.get(i)));
        }
        Double max = Collections.max(fies);
        Double min = Collections.min(fies);
        if (max * min >= 0) {
            result.add(P);
            return result;
        }
        for (int i = 0; i < n - 1; i++) {
            double fi = fies.get(i);
            double fi_1 = fies.get(i + 1);
            if ((Math.abs(fi) + Math.abs(fi_1)) == 0) {
                continue;
            }
            if (fi * fi_1 >= 0) {
                if (fi_1 > 0) {
                    Ll.add(P.get(i));
                    Ll.add(P.get(i + 1));
                } else {
                    Lr.add(P.get(i));
                    Lr.add(P.get(i + 1));
                }
            }
            Vector2D q = Vector2D.minus(P.get(i).multOnNumber(fi_1), P.get(i + 1).multOnNumber(fi)).divide(fi_1 - fi);
            L.add(Vector2D.multipl(Vector2D.minus(q, o), V));

            if (fi < 0) {
                Ll.add(P.get(i));
                Ll.add(q);
            }

            if (fi > 0) {
                Lr.add(P.get(i));
                Lr.add(q);
            }

            if (fi_1 > 0) {
                Lr.add(q);
                Lr.add(P.get(i + 1));
            }

            if (fi_1 < 0) {
                Ll.add(q);
                Ll.add(P.get(i + 1));
            }
        }

        Collections.sort(L);
        for (int k = 0; k < L.size() - 1; k++) {
            double Lk = L.get(k);
            double Lk_1 = L.get(k + 1);
            if (Lk_1 - Lk == 0) {
                continue;
            }
            Vector2D c = Vector2D.add(o, V.multOnNumber(Lk));
            Vector2D d = Vector2D.add(o, V.multOnNumber(Lk_1));
            double s = octTest(P, Vector2D.add(c, d).multOnNumber(0.5));
            if (s >= 0) {
                continue;
            }
            Lr.add(c);
            Lr.add(d);
            Ll.add(c);
            Ll.add(d);
        }
        List<List<Vector2D>> LLl = cont(Ll);
        List<List<Vector2D>> LLr = cont(Lr);
        result = union(LLl, LLr);
        return result;
    }

    private static List<List<Vector2D>> cont(List<Vector2D> Ls) {
        int m = 0;
        List<List<Vector2D>> LP = new ArrayList<>();
        while (!Ls.isEmpty()) {
            List<Vector2D> S = Ls;
            Vector2D a = Ls.get(0);
            Vector2D b = Ls.get(1);
            if(m == 0) {
                Line2D c = new Line2D.Double(a.x, a.y, b.x, b.y);
            }
        }
        return LP;
    }

    public static <T> List<T> union(List<T> list1, List<T> list2) {
        Set<T> set = new HashSet<>();

        set.addAll(list1);
        set.addAll(list2);

        return new ArrayList<>(set);
    }

    public static List<List<Vector2D>> triIteratePoly(List<Vector2D> P) throws Exception {
        List<List<Vector2D>> LP = new ArrayList<>();
        int d = dirTest(P);
        System.out.println(d);
        while (true) {
            System.out.println("WHILE OUTER");
            P = minPoly(P);
            int n = P.size() - 1;
            System.out.println(n + " ##############################################################################################");
            if (n == 3) {
                LP.add(P);
                return LP;
            }
            int v = 0;
            while (true) {

                System.out.println("WHILE INNER");
                List<Vector2D> T = new ArrayList<>();
                T.add(P.get(0));
                T.add(P.get(1));
                T.add(P.get(2));
                T.add(P.get(0));
                if (v++ == P.size() - 1) {
                    throw new Exception();
                }
                System.out.println(Arrays.deepToString(T.toArray()));
                if (d * nf2(P.get(0), P.get(1), P.get(2)) >= 0) {
                    LCShift(P);
                    continue;
                }
                boolean goToWhile = false;
                for (int i = 3; i < n; i++) {
                    int value = octTest(T, P.get(i));
                    System.out.println("VALUE = " + value);
                    if (value == 1) {
                        continue;
                    } else {
                        LCShift(P);
                        goToWhile = true;
                        break;
                    }
                }
                if (goToWhile) {
                    continue;
                }
                LP.add(T);
                P.remove(1);
                break;
            }

        }
    }

    public static int conv2Test(List<Vector2D> P, Vector2D q) {
        int r = 0, l = 0, e = 0;
        for (int i = 0; i < P.size() - 1; i++) {
            double f = nf2(P.get(i), P.get(i + 1), q);
            if (f < 0) {
                l = 1;
            }
            if (f > 0) {
                r = 1;
            }
            if (f == 0) {
                e = 1;
            }
            if (l * r != 0) {
                return 1;
            }
        }
        return e - 1;
    }

    public static int triPoly(List<Vector2D> P, List<List<Vector2D>> triangles) {
        System.out.println("####################################################");
        System.out.println("NEW VECTOR WITH SIZE" + " " + P.size());
        for (Vector2D asd : P) {
            System.out.println(asd.toString());
        }
        P = minPoly(P);
        System.out.println("MINIMAZE POLY" + " " + P.size());
        for (Vector2D asd : P) {
            System.out.println(asd.toString());
        }
        System.out.println("####################################################");
        int n = P.size() - 1;
        if (n == 3) {
            System.out.println("BECAUSE SIZE IS 4, WE ARE GOING TO END THIS STEP");
            System.out.println("NEW TRIANGLE ADDED ^^^^^^^^^");
            triangles.add(P);
            System.out.println("RETURN SIZE = " + triangles.size());
            return triangles.size();
        }

        int d = dirTest(P);
        System.out.println("DIR TEST RETURN " + d);
        int k = 3;
        System.out.println("CONV2(P) RETURN " + conv2(P));
        if (conv2(P) == 1) {
            System.out.println("BECAUSE POLYGON IS VIPUHLIY, GOING TO RETURN FUCTION");
            resurnFunction(P, k - 1, n, triangles);
            System.out.println("RETURN SIZE = " + triangles.size());
            return triangles.size();
        } else {
            System.out.println("POLYGON IS NOT VIPUHLIY");
            while (nf2(P.get(n - 1), P.get(0), P.get(1)) > 0) {
                System.out.println("LCShift");
                LCShift(P);
            }
            boolean flag = true;
            while (flag) {
                flag = false;
                while (true) {
                    if (k == n - 1) {
                        System.out.println("K = N - 1 GOING TO RETURN FUNC");
                        resurnFunction(P, k - 1, n, triangles);
                        System.out.println("RETURN SIZE = " + triangles.size());
                        return triangles.size();
                    }
                    double val = d * nf2(P.get(0), P.get(1), P.get(k - 1));
                    if (val >= 0) {
                        k++;
                    } else {
                        break;
                    }
                }
                for (int i = k + 1; i < n - 1; i++) {
                    double c = crossSegm(P.get(0), P.get(k - 1), P.get(i), P.get(i + 1));
                    if (c <= 0) {
                    } else {
                        System.out.println("CROSS SEGM RETURN POSITIVE VALUE, GOING TO RETURN FUNCTION");
                        k = i;
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    System.out.println("RETURN FUNCTION IN END GO TO");
                    resurnFunction(P, k - 1, n, triangles);
                    System.out.println("RETURN SIZE = " + triangles.size());
                    return triangles.size();
                }
            }
        }
        return 0;
    }

    public static List<Vector2D> cloneList(List<Vector2D> list) {
        List<Vector2D> clone = new ArrayList<>(list.size());
        for (Vector2D item : list) {
            clone.add(item.clone());
        }
        return clone;
    }

    public static double crossSegm(Vector2D a, Vector2D b, Vector2D c, Vector2D d) {
        System.out.println("CRESS SEGM IS OPEN FOR");
        System.out.println(a.toString());
        System.out.println(b.toString());
        System.out.println(c.toString());
        System.out.println(d.toString());
        double deltaXFirst = b.getX() - a.getX();
        double deltaYFirst = b.getY() - a.getY();

        double deltaXSec = c.getX() - d.getX();
        double deltaYSec = c.getY() - d.getY();

        double[][] matrix = new double[][]{{deltaXFirst, deltaYFirst}, {deltaXSec, deltaYSec}};
        int val = MatrixOperations.invert(matrix);
        if (val == 0) {
            return -1;
        }
        double[][] res = MatrixOperations.multiply(new double[][]{{c.getX() - a.getX(), c.getY() - a.getY()}}, matrix);
        System.out.println("CROSS SEGM RETURN " + (((0 <= res[0][0] && res[0][0] <= 1) ? 1 : 0) * ((0 <= res[0][1] && res[0][1] <= 1) ? 1 : 0)));
        return (((0 <= res[0][0] && res[0][0] <= 1) ? 1 : 0) * ((0 <= res[0][1] && res[0][1] <= 1) ? 1 : 0));
    }

    public static void resurnFunction(List<Vector2D> P, int k, int n, List<List<Vector2D>> triangles) {
        System.out.println("P size = " + P.size() + "k= " + k + "n= " + n);
        List<Vector2D> first = new ArrayList<>();
        first.add(P.get(0));
        //??
        for (int i = k; i <= n; i++) {
            first.add(P.get(i));
        }
        first.add(P.get(0));
        triPoly(first, triangles);
        first.clear();
        first.add(P.get(0));
        for (int i = 1; i <= k; i++) {
            first.add(P.get(i));
        }
        first.add(P.get(0));
        triPoly(first, triangles);
    }

    public static void LCShift(List<Vector2D> P) {
//        System.out.println("LCShift IS OPEN FOR ");
//        Arrays.deepToString(P.toArray());
        P.add(P.get(1));
        P.remove(0);
//        System.out.println("END OF LCShift");
//        Arrays.deepToString(P.toArray());
    }

    public static List<Vector2D> minPoly(List<Vector2D> P) {
//        System.out.println("MIN POLY IS OPEN");
        int m = 0;
        for (int i = 1; i < P.size(); i++) {
            Vector2D V = Vector2D.minus(P.get(i), P.get(m));
            if (Vector2D.multipl(V, V) == 0) {
                continue;
            }
            if ((m == 0) || (nf2(P.get(m - 1), P.get(m), P.get(i)) != 0)) {
                m++;
                P.get(m).setVector2D(P.get(i));
                continue;
            }
            Vector2D W = Vector2D.minus(P.get(m), P.get(m - 1));
            if (Vector2D.multipl(V, W) > 0) {
                P.get(m).setVector2D(P.get(i));
                continue;
            }
            if (Vector2D.multipl(Vector2D.minus(P.get(i), P.get(m - 1)), W) > 0) {
                continue;
            }
            P.get(m - 1).setVector2D(P.get(m));
            P.get(m).setVector2D(P.get(i));
        }
        List<Vector2D> res = new ArrayList<>();
        if (nf2(P.get(m - 1), P.get(0), P.get(1)) == 0) {
            res.add(P.get(1));
            for (int j = 2; j < m; j++) {
                res.add(P.get(j));
            }
            res.add(P.get(1));
        } else {
            for (int j = 0; j <= m; j++) {
                res.add(P.get(j));
            }
        }
//        System.out.println("MIN POLY RETURN");
        return res;
    }

    public static int conv2(List<Vector2D> P) {
        int n = P.size() - 1;
        double s = nf2(P.get(P.size() - 2), P.get(0), P.get(1));
        double f;
        //////????
        for (int i = 1; i < n; i++) {
            f = nf2(P.get(i - 1), P.get(i), P.get(i + 1));
            if (s * f > 0) {
                continue;
            }
            if (s * f == 0) {
                return -1;
            }
            return 0;
        }

        return 1;
    }

    public static int dirTest(List<Vector2D> P) {
//        System.out.println("DIR TEST IS OPEN");
//        System.out.println("P0 = " + P.get(0).toString() + "P(P.size -1) " + P.get(P.size() - 1).toString());
        Vector2D Qgr = Vector2D.add(P.get(0), P.get(P.size() - 1));
//        System.out.println("Qgr = " + Qgr.toString());
        double sum = 0;
        for (int i = 0; i < P.size() - 1; i++) {
            sum += ang(Vector2D.minus(P.get(i), Qgr), Vector2D.minus(P.get(i + 1), Qgr));
        }
//        System.out.println("SUM IN DIR TEST = " + sum);
        return sign(sum);
    }

    public static double ang(Vector2D V, Vector2D W) {
        boolean positive = false;
        if (MatrixOperations.det(new double[][]{{V.x, V.y}, {W.x, W.y}}) >= 0) {
            positive = true;
        }
        return (positive ? 1 : -1) * Math.acos(Vector2D.multipl(V, W) / (V.module() * W.module()));
    }

    public static double nfV(Vector2D o, Vector2D V, Vector2D p) {
        double[][] first = new double[1][2];
        first[0][0] = p.getX() - o.getX();
        first[0][1] = p.getY() - o.getY();

        double[][] res = new double[2][2];
        res[0][0] = first[0][0];
        res[0][1] = first[0][1];
        res[1][0] = V.x;
        res[1][1] = V.y;
        double val = MatrixOperations.det(res);
        return val;
    }

    public static double nf2(Vector2D a, Vector2D b, Vector2D p) {
//        System.out.println("NF2 IS OPEN  FOR");
//        System.out.println(a.toString());
//        System.out.println(b.toString());
//        System.out.println(p.toString());
        double[][] first = new double[1][2];
        first[0][0] = p.getX() - a.getX();
        first[0][1] = p.getY() - a.getY();

        double[][] second = new double[2][1];
        second[0][0] = b.getY() - a.getY();
        second[1][0] = -(b.getX() - a.getX());
        double[][] res = new double[2][2];
        res[0][0] = first[0][0];
        res[0][1] = first[0][1];
        res[1][0] = b.x - a.x;
        res[1][1] = b.y - a.y;
        double val = MatrixOperations.det(res);
//        System.out.println("NF2 RETURN " + val);
        return val;
//        double[][] asd = MatrixOperations.multiply(first, second);
////        System.out.println("NF2 RETURN " + asd[0][0]);
//        return asd[0][0];
    }

    public static int sign(double x) {
        System.out.println("SIGN IS OPEN FOR " + x);
        if (x > 0) {
            return 1;
        }
        if (x < 0) {
            return -1;
        }
        return 0;
    }

    public static int octTest(List<Vector2D> P, Vector2D q) {
        double s = 0;
        int w = 0;
        for (int i = 0; i < P.size(); i++) {
            Vector2D pp = P.get(i);
            int v = oct(pp.getX() - q.getX(), pp.getY() - q.getY());
            if (v == 0) {
                return 0;
            }
            if (i == 0) {
                w = v;
                continue;
            }
            int delta = v - w;
            if (Math.abs(delta) < 4) {
                s += delta;
                w = v;
                continue;
            }

            if (Math.abs(delta) > 4) {
                delta = delta - 8 * sign(delta);
                s += delta;
                w = v;
                continue;
            }
            double f = nf2(P.get(i - 1), pp, q);
            if (f == 0) {
                return 0;
            }
            delta = -4 * sign(f);
            s += delta;
            w = v;
        }
        return 1 - 2 * sign((s));
    }

    private static int oct(double x, double y) {
        if (x == 0 && y == 0) {
            return 0;
        }
        if (0 <= y && y < x) {
            return 1;
        }
        if (0 < x && x <= y) {
            return 2;
        }
        if (-x <= y && y < 0) {
            return 8;
        }
        if (0 <= x && x < -y) {
            return 7;
        }
        if (0 < y && y <= -x) {
            return 4;
        }
        if (-y < x && x <= 0) {
            return 3;
        }
        if (x < y && y <= 0) {
            return 5;
        }
        if (y <= x && x < 0) {
            return 6;
        }

        return -100;
    }
}
