/*
 * PROJECT III: TriMatrix.java
 *
 * This file contains a template for the class TriMatrix. Not all methods are
 * implemented. Make sure you have carefully read the project formulation
 * before starting to work on this file. You will also need to have completed
 * the Matrix class.
 *
 * Remember not to change the names, parameters or return types of any
 * variables in this file!
 *
 * The function of the methods and instance variables are outlined in the
 * comments directly above them.
 */

public class TriMatrix extends Matrix {
    /**
     * An array holding the diagonal elements of the matrix.
     */
    private double[] diag;

    /**
     * An array holding the upper-diagonal elements of the matrix.
     */
    private double[] upper;

    /**
     * An array holding the lower-diagonal elements of the matrix.
     */
    private double[] lower;
    
    /**
     * Constructor function: should initialise m and n through the Matrix
     * constructor and set up the data array.
     *
     * @param N  The dimension of the array.
     */
    public TriMatrix(int N) {
	super(N,N);
	if(N<1){
	  throw new MatrixException("Matrix cannot have dimension less than or equal to 0");
	}
	diag = new double[N];
	upper = new double[N-1];
	lower = new double[N-1];
    }
    
    /**
     * Getter function: return the (i,j)'th entry of the matrix.
     *
     * @param i  The location in the first co-ordinate.
     * @param j  The location in the second co-ordinate.
     * @return   The (i,j)'th entry of the matrix.
     */
    public double getIJ(int i, int j) {
	if(i>=0 && i<m && j>=0 && j<n){
	  if(i==j){
	    return diag[i];
	  }else if(j==i+1){
	    return upper[i];
	  }else if(i==j+1){
	    return lower[j];
	  }else{
	    return 0;
	  }
	}else{
	  throw new MatrixException("Cannot return element as since it is not in the matrix");
	}
    }
    
    /**
     * Setter function: set the (i,j)'th entry of the data array.
     *
     * @param i    The location in the first co-ordinate.
     * @param j    The location in the second co-ordinate.
     * @param val  The value to set the (i,j)'th entry to.
     */
    public void setIJ(int i, int j, double val) {
	if(i>=0 && i<m && j>=0 && j<n){
	  if(i==j){
	    diag[i] = val;
	  }else if(j==i+1){
	    upper[i] = val;
	  }else if(i==j+1){
	    lower[j] = val;
	  }else{
	    throw new MatrixException("Cannot set non tridiagonal element.");
	  }      
	}else{
	  throw new MatrixException("Cannot set element as it is not in the matrix");	
	}
    }
    
    /**
     * Return the determinant of this matrix.
     *
     * @return The determinant of the matrix.
     */
    public double determinant() {
	double det = 1;
	TriMatrix LU = this.decomp();
	for(int i=0;i<this.n;i++){
	    det = det*LU.getIJ(i,i);
	}
	return det;
    }
    
    /**
     * Returns the LU decomposition of this matrix. See the formulation for a
     * more detailed description.
     * 
     * @return The LU decomposition of this matrix.
     */
    public TriMatrix decomp() {
        TriMatrix LU = new TriMatrix(this.n);
	LU.setIJ(0,0,this.getIJ(0,0));
	for(int u=0;u<this.n-1;u++){
	    LU.setIJ(u,u+1,this.getIJ(u,u+1));
	}
	for(int d=0;d<this.n-1;d++){
	    LU.setIJ(d+1,d,this.getIJ(d+1,d)/LU.getIJ(d,d)); 
	    LU.setIJ(d+1,d+1,this.getIJ(d+1,d+1)-(LU.getIJ(d+1,d)*LU.getIJ(d,d+1)));
	}
	return LU;
    }

    /**
     * Add the matrix to another matrix A.
     *
     * @param A  The Matrix to add to this matrix.
     * @return   The sum of this matrix with the matrix A.
     */
    public Matrix add(Matrix A){
	if(A instanceof GeneralMatrix){
	  return A.add(this);
	}else if(A instanceof TriMatrix){
	  double total=0;
	  if(this.n==A.n){
	    TriMatrix rMatrix = new TriMatrix(this.n);
	    for(int u=0;u<this.n-1;u++){
		total=this.getIJ(u,u+1)+A.getIJ(u,u+1);
		rMatrix.setIJ(u,u+1,total);
	    }
	    for(int d=0;d<this.n;d++){
		total=this.getIJ(d,d)+A.getIJ(d,d);
		rMatrix.setIJ(d,d,total);
	    }
	    for(int l=0;l<this.n-1;l++){
		total=this.getIJ(l+1,l)+A.getIJ(l+1,l);
		rMatrix.setIJ(l+1,l,total);
	    }			
	    return rMatrix;
	  }else{
		throw new MatrixException("Dimensions are different, so cannot add");
	  }	
	}else{
	      throw new MatrixException("Unknown Matrix Type");
	}
    }
    
    /**
     * Multiply the matrix by another matrix A. This is a _left_ product,
     * i.e. if this matrix is called B then it calculates the product BA.
     *
     * @param A  The Matrix to multiply by.
     * @return   The product of this matrix with the matrix A.
     */
    public Matrix multiply(Matrix A) {
	double rowt;
	if(this.n==A.m){
	  GeneralMatrix rMatrix = new GeneralMatrix(this.m,A.n);
	  for(int r=0;r<this.m;r++){
	      for(int c=0;c<A.n;c++){
		  rowt=0;
		  for(int a=0;a<this.n;a++){
		      rowt=rowt+(this.getIJ(r,a)*A.getIJ(a,c));
		  }	
		  rMatrix.setIJ(r,c,rowt);
	      }
	  }
	  return rMatrix;
	}else{
	      throw new MatrixException("Cannot Multiply, as width of first isnt equal to height of second");
	}
    }
    
    /**
     * Multiply the matrix by a scalar.
     *
     * @param a  The scalar to multiply the matrix by.
     * @return   The product of this matrix with the scalar a.
     */
    public Matrix multiply(double a) {
	TriMatrix rMatrix = new TriMatrix(this.n);
	for(int u=0;u<this.n-1;u++){
	    rMatrix.setIJ(u,u+1,this.getIJ(u,u+1)*a);
	}
	for(int d=0;d<this.n;d++){
	    rMatrix.setIJ(d,d,this.getIJ(d,d)*a);
	}
	for(int l=0;l<this.n-1;l++){
	    rMatrix.setIJ(l+1,l,this.getIJ(l+1,l)*a);
	}
	return rMatrix;
    }

    /**
     * Populates the matrix with random numbers which are uniformly
     * distributed between 0 and 1.
     */
    public void random() {
	for(int u=0;u<this.n-1;u++){
	    this.setIJ(u,u+1,Math.random());
	}
	for(int d=0;d<this.n;d++){
	    this.setIJ(d,d,Math.random());
	}
	for(int l=0;l<this.n-1;l++){
	    this.setIJ(l+1,l,Math.random());
	}
    }
    
    /*
     * Your tester function should go here.
     */
    public static void main(String[] args) {
	Matrix A = new TriMatrix(3);
	A.setIJ(0,0,1);
	A.setIJ(0,1,2);
	A.setIJ(1,0,3);
	A.setIJ(1,1,1);
	A.setIJ(1,2,1);
	A.setIJ(2,1,2);
	A.setIJ(2,2,3);
	System.out.println("(0,0) element of A is: "+A.getIJ(0,0));        
	System.out.println("TriMatrix A: \n"+A.toString());
	System.out.println("Determinant of A is: "+A.determinant()+"\n");

	Matrix B = new TriMatrix(3);
	B.setIJ(0,0,4);
	B.setIJ(0,1,1);
	B.setIJ(1,0,2);
	B.setIJ(1,1,2);
	B.setIJ(1,2,2);
	B.setIJ(2,1,4);
	B.setIJ(2,2,1);
	System.out.println("TriMatrix B: \n"+B.toString());
        System.out.println("Determinant of B is: "+B.determinant()+"\n");
	
	Matrix C = A.add(B);
	System.out.println("A+B: \n"+C.toString());

	Matrix D = new GeneralMatrix(3,3);
	D.setIJ(0,0,1);
	D.setIJ(0,1,1);
	D.setIJ(1,0,1);
	D.setIJ(2,0,1);
	D.setIJ(0,2,1);
	D.setIJ(1,1,1);
	D.setIJ(1,2,1);
	D.setIJ(2,1,1);
	D.setIJ(2,2,1);	
	System.out.println("GeneralMatrix D: \n"+D.toString());
		
	Matrix AD = A.add(D);
	System.out.println("A+D: \n"+AD.toString());

	Matrix E = A.multiply(3);
	System.out.println("3*A: \n"+E.toString());
	
	Matrix F = A.multiply(B);
	System.out.println("A*B: \n"+F.toString());

	Matrix G = A.multiply(D);
	System.out.println("A*D: \n"+G.toString());

	Matrix H = new TriMatrix(4);
	H.random();
	System.out.println("Random TriMatrix: \n"+H.toString());

	System.out.println("Trying to premultiply 3x3 trimatrix with 5X5 trimatrix: ");
	try{
	    Matrix a = new TriMatrix(5);
	    a.random();
	    Matrix b = new TriMatrix(3);
	    b.random();
	    Matrix c = b.multiply(a);
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}
	
	System.out.println("\nTrying to add 11x11 trimatrix to 4X4 trimatrix: ");
	try{
	    Matrix a1 = new TriMatrix(11);
	    a1.random();
	    Matrix b1 = new TriMatrix(4);
	    b1.random();
	    Matrix c1 = b1.add(a1);
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}				
	
	System.out.println("\nTrying to set a non upper/lower/centre diagonal element: ");
	try{
	    Matrix a2 = new TriMatrix(10);
	    a2.setIJ(10,0,2.5);
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}
	
	System.out.println("\nTrying to Create Matrix with one dimension zero ");
	try{
	    Matrix a3 = new TriMatrix(0);	
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}

    }
}
