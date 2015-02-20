package ru.saintcat.client.interpolcurves;

import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.ResourceBundle;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.geometry.Point2D;
import javafx.scene.control.Label;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.AnchorPane;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Line;
import javafx.scene.shape.Polygon;

public class CurveController implements Initializable {

    @FXML
    private AnchorPane paintPanel;
    @FXML
    private Label coordinatsLabel;
    @FXML
    private Label resultLabel;
    private ObservableList<Circle> circles = FXCollections.observableArrayList();
    private List<Line> lines = new ArrayList<>();
    private List<Vector2D> points = new ArrayList<>();
    private List<Polygon> pols = new ArrayList<>();
    private Circle point;
    private static final int CIRCLE_SIZE = 6;

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        newPolAction(null);
        paintPanel.addEventFilter(MouseEvent.MOUSE_MOVED, new EventHandler<MouseEvent>() {

            @Override
            public void handle(MouseEvent t) {
                coordinatsLabel.setText(t.getX()+ " " + t.getY());
            }
        });
    }
    private double shift = 0.2;

    private void triPolyTick() {
        List<List<Vector2D>> result = new ArrayList<>();
        //        TrianglePoly.triPoly(test, result);
//        points.clear();
//        points.add(new Vector2D(100, 100));
//        points.add(new Vector2D(100, 300));
//        points.add(new Vector2D(300, 300));
//        points.add(new Vector2D(300, 100));
//        points.add(new Vector2D(100, 100));
        result = TrianglePoly.cutLine(points, new Vector2D(200, 50), new Vector2D(500, 400));
        resultLabel.setText("Результат " + result.size() + " треугольника(ов)");
        double shift = 0;
        int f = 0;
        Random rand = new Random();
        paintPanel.getChildren().removeAll(pols);
        pols.clear();
        Line p = new Line(200, 50, 200 + 500, 50 + 400);
        paintPanel.getChildren().add(p);
        for (List<Vector2D> r : result) {
            int rCol1 = rand.nextInt(256);
            int rCol2 = rand.nextInt(256);
            int rCol3 = rand.nextInt(256);
            Polygon polygon = new Polygon();
            polygon.getPoints().addAll(getPoints(r));
            polygon.setFill(Color.rgb(rCol1, rCol2, rCol3));
            pols.add(polygon);
            paintPanel.getChildren().add(polygon);
            System.out.println("");
            f++;
        }
    }

    private int octTest(List<Vector2D> P, Point2D q) {
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

    private double nf2(Vector2D a, Vector2D b, Point2D p) {
        double[][] first = new double[1][2];
        first[0][0] = p.getX() - a.getX();
        first[0][1] = p.getY() - a.getY();

        double[][] second = new double[2][1];
        second[0][0] = b.getY() - a.getY();
        second[1][0] = -(b.getX() - a.getX());

        double[][] res = MatrixOperations.multiply(first, second);
        System.out.println(res.length + " " + res[0].length);
        return res[0][0];
    }

    private int sign(double x) {
        if (x > 0) {
            return 1;
        }
        if (x < 0) {
            return -1;
        }
        return 0;
    }

    private int oct(double x, double y) {
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

    private List<Vector2D> apoly(Vector2D T, int a, int b, int O) {
        List<Vector2D> res = new ArrayList<>();
        double fi = 0;
        while (fi < 360) {
            double d = a + Rand.nextInt(b - a);
            double dx = T.x + d * Math.cos(fi * Math.PI / 180);
            double dy = T.y + d * Math.sin(fi * Math.PI / 180);
            System.out.println(dx + " " + dy);
            res.add(new Vector2D(dx, dy));
            System.out.println(fi);
            fi = fi + Rand.nextInt(O);
        }
        return res;
    }
    private static Random Rand = new Random();

    @FXML
    void newPolAction(ActionEvent event) {
        points.clear();
        paintPanel.getChildren().removeAll(lines);
        lines.clear();
        List<Vector2D> pol = apoly(new Vector2D(paintPanel.getPrefWidth() / 2 - 160, paintPanel.getPrefHeight() / 2), 10, 200, 20);
        points.addAll(pol);
        points.add(pol.get(0));
        for (int i = 0; i < points.size() - 1; i++) {
            Line line = new Line(points.get(i).getX(), points.get(i).getY(), points.get(i + 1).getX(), points.get(i + 1).getY());
            lines.add(line);
            paintPanel.getChildren().add(line);
        }
        triPolyTick();
    }

    public static List<Double> getPoints(List<Vector2D> vect) {
        List<Double> res = new ArrayList<>();
        int z = 0;
        for (Vector2D v : vect) {
            res.add(v.x + 350);
            res.add(v.y);
        }
        return res;
    }
}
