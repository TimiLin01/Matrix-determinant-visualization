/*
 * PROJECT III: Project3.java
 *
 * This file contains a template for the class Project3. None of methods are
 * implemented. Make sure you have carefully read the project formulation
 * before starting to work on this file. You will also need to have completed
 * the Matrix class, as well as GeneralMatrix and TriMatrix.
 *
 * Remember not to change the names, parameters or return types of any
 * variables in this file!
 *
 * The function of the methods and instance variables are outlined in the
 * comments directly above them.
 */

public class Project3 {
    /**
     * Calculates the variance of the distribution defined by the determinant
     * of a random matrix. See the formulation for a detailed description.
     *
     * @param m           The matrix object that will be filled with random
     *                    samples.
     * @param numSamples  The number of samples to generate when calculating 
     *                    the variance. 
     * @return            The variance of the distribution.
     */
    public static double matVariance(Matrix m, int numSamples) {
	double det;
	double sumdet = 0;
	double sumdetsq = 0;
	for(int i=0;i<numSamples;i++){
	    m.random();
	    det=m.determinant();
	    sumdet = sumdet+det;
	    sumdetsq=sumdetsq+Math.pow(det,2);
	}
	return (sumdetsq/numSamples-Math.pow((sumdet/numSamples),2));
    }
    
    /**
     * This function should calculate the variances of matrices for matrices
     * of size 2 <= n <= 50. See the formulation for more detail.
     */
    public static void main(String[] args) {
	GeneralMatrix gmat;
	TriMatrix tmat;
	for(int n=2;n<=50;n++){
	    gmat = new GeneralMatrix(n,n);
	    tmat = new TriMatrix(n);
	    System.out.println(n+" "+matVariance(gmat,15000)+" "+matVariance(tmat,150000));
	}
    }
}
