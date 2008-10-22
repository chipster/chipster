/*
 * Matrix.java
 *
 * Created on 24. elokuuta 2006, 0:40
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package fi.csc.microarray.client.visualisation.methods.threed;

/**
 * From the Viski project (http://www.cs.helsinki.fi/group/viski/).
 *
 * @author esa
 */
public class Matrix {
    
    public double array[][];
    
    /**
     * Creates a new instance of Matrix
     * @param m 
     * @param n 
     */
    public Matrix(int m, int n) {
        this.array = new double[m][n];
    }
    
    /**
     * 
     * @param m 
     * @param n 
     * @param v 
     */
    public Matrix(int m, int n, double v) {
        this(m, n);
        for (int i=0; i < array.length; ++i)
            java.util.Arrays.fill(array[i], v);
    }
    
    /**
     * 
     * @param m 
     * @param n 
     * @param v 
     */
    public void set(int m, int n, double v) {
        this.array[m][n] = v;
    }
    
    /**
     * 
     * @param m 
     * @param n 
     * @return 
     */
    public double get(int m, int n) {
        return this.array[m][n];
    }
    
    /**
     * 
     * @param that 
     * @return 
     */
    public Matrix times(Matrix that) {
        if (this.array[0].length != that.array.length)
            throw new IllegalArgumentException("Matrix dimension do not allow multiplication.");
        
        int m = this.array.length;
        int n = this.array[0].length;
        int a = that.array[0].length;
        double[] tmp = new double[n];
                
        Matrix r = new Matrix(m, a);
        
        for (int j=0; j < a; ++j) {
            for (int i=0; i < n; ++i)
                tmp[i] = that.array[i][j];
            for (int i=0; i < m; ++i) {
                r.array[i][j] = times(this.array[i], tmp);
            }
        }
        return r;
    }
    
    private double times(double[] a, double[] b) {
        if (a.length != b.length)
            throw new IllegalArgumentException("Vector sizes are different.");
        
        double sum = 0;
        for (int i=0; i < a.length; ++i) {
            sum += a[i] * b[i];
        }
        
        return sum;
    }
}
