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
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import org.apache.commons.collections4.CollectionUtils;

/**
 *
 * @author Roman Chernyshev
 */
public class TrianglePoly {

    public static List<List<Vector2D>> cutLine(List<Vector2D> P, Vector2D o, Vector2D V) {
        int n = P.size();
        System.out.println(V);
        V = V.normal();
         System.out.println(V);
        List<List<Vector2D>> result = new ArrayList<>();
        List<Double> L = new ArrayList<>();
        List<MyLine2D> Ll = new ArrayList<>();
        List<MyLine2D> Lr = new ArrayList<>();
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
                    Lr.add(new MyLine2D(P.get(i), P.get(i + 1)));
                } else {
                    Ll.add(new MyLine2D(P.get(i), P.get(i + 1)));
                }
                continue;
            }
            Vector2D q = Vector2D.minus(P.get(i).multOnNumber(fi_1), P.get(i + 1).multOnNumber(fi)).divide(fi_1 - fi);
            L.add(Vector2D.multipl(Vector2D.minus(q, o), V));
            if (fi < 0) {
                Ll.add(new MyLine2D(P.get(i), q));
            }
            if (fi > 0) {
                Lr.add(new MyLine2D(P.get(i), q));
            }
            if (fi_1 > 0) {
                Lr.add(new MyLine2D(q, P.get(i + 1)));
            }
            if (fi_1 < 0) {
                Ll.add(new MyLine2D(q, P.get(i + 1)));
            }
        }
        System.out.println("##########" + Ll);
        System.out.println("##########" + Lr);
        System.out.println("##########" + L);
        Collections.sort(L);
        System.out.println("##########" + L);
        for (int k = 0; k < L.size() - 1; k++) {
            double Lk = L.get(k);
            System.out.println("Lk " + Lk);
            double Lk_1 = L.get(k + 1);
            System.out.println("Lk + 1 " + Lk_1);
            if (Lk_1 - Lk == 0) {
                System.out.println("CONTINUE1");
                continue;
            }
            Vector2D c = Vector2D.add(o, V.multOnNumber(Lk));
            System.out.println("c " + c);
            Vector2D d = Vector2D.add(o, V.multOnNumber(Lk_1));
            System.out.println("d " + d);
            System.out.println("middle point " + Vector2D.add(c, d).multOnNumber(0.5));
            double s = octTest(P, Vector2D.add(c, d).multOnNumber(0.5));
            System.out.println("oct test = " + s);
            if (s > 0) {
                System.out.println("CONTINUE2");
                continue;
            }
            Lr.add(new MyLine2D(c, d));
            Ll.add(new MyLine2D(c, d));
        }
        System.out.println("BEFORE CONT");
        System.out.println("##########" + Ll);
        System.out.println("##########" + Lr);
        List<List<Vector2D>> LLl = cont(transform(Ll));
        System.out.println("ZZZZ##########" + LLl);
        List<List<Vector2D>> LLr = cont(transform(Lr));
        System.out.println("ZZZZ##########" + LLr);
        List<List<Vector2D>> ff = union(LLl, LLr);
        result.addAll(ff);
        return result;
    }

    private static List<Vector2D> transform(List<MyLine2D> l) {
        List<Vector2D> res = new ArrayList<>();
        for (MyLine2D ln : l) {
            res.add(ln.getFirstPoint());
            res.add(ln.getSecondPoint());
        }
        return res;
    }

    private static List<List<Vector2D>> cont(List<Vector2D> Ls) {
        System.out.println(Ls);
//        int m = 0;
        List<Vector2D> C = new ArrayList<>();
        List<List<Vector2D>> LP = new ArrayList<>();
        while (!Ls.isEmpty()) {
//            System.out.println("LS NOT EMPTY");
            List<Vector2D> S = new ArrayList<>(Ls.size());
            CollectionUtils.addAll(S, Ls);
//            System.out.println("S AFTER COPY = " + S);
            int f = 1;
            Vector2D a = null;
            Vector2D b = null;
            while (f == 1) {
//                System.out.println("F == 1");
                a = Ls.get(0);
                b = Ls.get(1);
                if (C.size() == 0) {
                    C.add(a);
                    C.add(b);
                    break;
                }
                if (Vector2D.minus(C.get(C.size() - 1), a).isZero()) {
                    System.out.println("ZZZZZZZZZZZZz");
                    C.add(b);
                    break;
                }
                if (Vector2D.minus(C.get(C.size() - 1), b).isZero()) {
                    System.out.println("FFFFFFFFFFFFFF");
                    C.add(a);
                    break;
                }
                LCShift(Ls);
//                System.out.println("LS = " + Ls.size());
//                System.out.println("S = " + S.size());
                if (collectEquals(Ls, S)) {
//                    System.out.print("F = 0");

                    C.remove(C.size() - 1);
//                    m = m - 1 + (f = 0);
                    break;
                }
            }
            if (f == 1) {
                Ls.remove(a);
                Ls.remove(b);
//                m = C.size() - 1;
                if (C.size() >= 3) {
                    for (int j = 0; j <= C.size() - 3; j++) {
//                        System.out.println("INNER");
//                        System.out.println("Cm = " + C.get(m));
//                        System.out.println("Cj = " + C.get(j));
                        if (Vector2D.minus(C.get(C.size() - 1), C.get(j)).isZero()) {
                            if (j != 0) {
                                for (int k = 1; k <= j; k++) {
                                    Ls.add(C.get(0));
                                    Ls.add(C.get(1));
                                    C.remove(0);
                                }
                            }
                            System.out.println("ADDDDDDDDDDDDDDD");
                            LP.add(C);
//                            m = 0;
                            break;
                        }
                    }

                }
            }
        }
        return LP;
    }

    private static boolean collectEquals(List<Vector2D> f, List<Vector2D> s) {
        for (int i = 0; i < f.size(); i++) {
            if (f.get(i).equals(s.get(i))) {
            } else {
                System.out.println("RETURN FALSE");
                return false;
            }
        }
        System.out.println("RETURN TRUE");
        return true;
    }

    public static <T> List<T> union(List<T> list1, List<T> list2) {
        Set<T> set = new HashSet<>();

        set.addAll(list1);
        set.addAll(list2);

        return new ArrayList<>(set);
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

    public static void LCShift(List<Vector2D> P) {
//        System.out.println("LCShift IS OPEN FOR ");
//        Arrays.deepToString(P.toArray());
        P.add(P.get(0));
        P.add(P.get(1));
        P.remove(0);
        P.remove(0);
//        System.out.println("END OF LCShift");
//        Arrays.deepToString(P.toArray());
    }

    public static void LCShift2(List<Vector2D> P) {
//        System.out.println("LCShift IS OPEN FOR ");
//        Arrays.deepToString(P.toArray());
        P.add(P.get(2));
        P.add(P.get(3));
        P.remove(0);
        P.remove(1);
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

    private static double nf2(Vector2D a, Vector2D b, Vector2D p) {
        double[][] first = new double[1][2];
        first[0][0] = p.getX() - a.getX();
        first[0][1] = p.getY() - a.getY();

        double[][] second = new double[2][1];
        second[0][0] = b.getY() - a.getY();
        second[1][0] = -(b.getX() - a.getX());

        double[][] res = MatrixOperations.multiply(first, second);
        return res[0][0];
    }

    public static int sign(double x) {
        if (x > 0) {
            return 1;
        }
        if (x < 0) {
            return -1;
        }
        return 0;
    }

    private static int octTest(List<Vector2D> P, Vector2D q) {
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
        System.out.print(1 - 2 * sign((s)));
        return 1 - 2 * sign((Math.abs(s)));
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