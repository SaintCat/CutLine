/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ru.saintcat.client.interpolcurves;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

/**
 *
 * @author ROMAN
 */
public class MyLine2D extends Line2D.Double {

    public MyLine2D(Vector2D v, Vector2D d) {
        super(v.x, v.y, d.x, d.y);
    }
    
    
    public Vector2D getFirstPoint() {
        return new Vector2D(x1,y1);
    }
    
    public Vector2D getSecondPoint() {
        return new Vector2D(x2,y2);
    }

    @Override
    public String toString() {
        return "MyLine2D{" + this.x1 + " " + y1 + " " + x2 + " " + y2 + '}';
    }
}
