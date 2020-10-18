import java.util.Arrays;

public class GeneralMatrix extends Matrix {
    /* This instance variable stores the elements of the matrix.*/
    private double[][] data;

    /* Constructor function: should initialise m and n through the Matrix
     * constructor and set up the data array.
     * @param m  The first dimension of the array.
     * @param n  The second dimension of the array.*/
    public GeneralMatrix(int m, int n) throws MatrixException {
   	super(m,n);
	if(m<1 || n<1){
	   throw new MatrixException("Matrix cannot have dimension less than or equal to 0");
	}
	data = new double[m][n];     
    }

    /* Constructor function. This is a copy constructor; it should create a
     * copy of the matrix A.
     * @param A  The matrix to create a copy of.*/
    public GeneralMatrix(GeneralMatrix A) {
	super(A.m,A.n);
	data = new double[A.m][A.n];
	for(int i=0;i<A.m;i++){
	    for(int j=0;j<A.n;j++){
		data[i][j] = A.getIJ(i,j);
	    }
	}
    }
    
    /* Getter function: return the (i,j)'th entry of the matrix.
     * @param i  The location in the first co-ordinate.
     * @param j  The location in the second co-ordinate.
     * @return   The (i,j)'th entry of the matrix.*/
    public double getIJ(int i, int j) {
        if(i>=0 && i<m && j>=0 && j<n){
	   return data[i][j];
	}else{
	   throw new MatrixException("Cannot return element since it is not in the matrix");
	}
    }
    
    /* Setter function: set the (i,j)'th entry of the data array.
     * @param i    The location in the first co-ordinate.
     * @param j    The location in the second co-ordinate.
     * @param val  The value to set the (i,j)'th entry to.*/
    public void setIJ(int i, int j, double val) {
        if(i>=0 && i<m && j>=0 && j<n){
            data[i][j]=val;
        }else{
            throw new MatrixException("cannot return element since it is not in the matrix");
        }
    }
    
    /* Return the determinant of this matrix.
     * @return The determinant of the matrix.*/
    public double determinant() {
	double[] d = new double[1];
	double det = 1;
	GeneralMatrix LUCombined = this.decomp(d);
	for(int i=0;i<this.n;i++){
	    det = det*LUCombined.getIJ(i,i);
	}
	det=det*d[0];
	return det;
    }

    /* Add the matrix to another matrix A.
     * @param A  The Matrix to add to this matrix.
     * @return   The sum of this matrix with the matrix A.*/
    public Matrix add(Matrix A) {
	if(this.n==A.n && this.m==A.m){
	    GeneralMatrix rMatrix = new GeneralMatrix(this.m,this.n);
	    for(int i=0;i<this.m;i++){
		for(int j=0;j<this.n;j++){
		    rMatrix.setIJ(i,j,(A.getIJ(i,j)+this.getIJ(i,j)));
		}
	    }
	return rMatrix;
	}else{
	 throw new MatrixException("Dimensions are different, so cannot add");
	}
    }
    
    /* Multiply the matrix by another matrix A. This is a _left_ product,
     * i.e. if this matrix is called B then it calculates the product BA.
     * @param A  The Matrix to multiply by.
     * @return   The product of this matrix with the matrix A.*/
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
		throw new MatrixException("Cannot Multiply, as width of first is not equal to height of second");
	  }
    }

    /* Multiply the matrix by a scalar.
     * @param a  The scalar to multiply the matrix by.
     * @return   The product of this matrix with the scalar a.*/
    public Matrix multiply(double a) {
	GeneralMatrix rMatrix = new GeneralMatrix(this.m,this.n); 
	for(int i=0;i<this.m;i++){
	    for(int j=0;j<this.n;j++){
		rMatrix.setIJ(i,j,getIJ(i,j)*a);
	    }
	}
	return rMatrix;
    }


    /* Populates the matrix with random numbers which are uniformly
     * distributed between 0 and 1.*/
    public void random() {
	for(int i=0;i<this.m;i++){
	    for(int j=0;j<this.n;j++){
		this.setIJ(i,j,Math.random());
	    }
	}
    }

    /* Returns the LU decomposition of this matrix; i.e. two matrices L and U
     * so that A = LU, where L is lower-diagonal and U is upper-diagonal. 
     * @param d  An array of length 1. On exit, the value contained in here
     *           will either be 1 or -1, which you can use to calculate the
     *           correct sign on the determinant.
     * @return   The LU decomposition of the matrix.*/
    public GeneralMatrix decomp(double[] d) {
        if (n != m)
            throw new MatrixException("Matrix is not square");
        if (d.length != 1)
            throw new MatrixException("d should be of length 1");
        
        int           i, imax = -10, j, k; 
        double        big, dum, sum, temp;
        double[]      vv   = new double[n];
        GeneralMatrix a    = new GeneralMatrix(this);
        
        d[0] = 1.0;
        
        for (i = 1; i <= n; i++) {
            big = 0.0;
            for (j = 1; j <= n; j++)
                if ((temp = Math.abs(a.data[i-1][j-1])) > big)
                    big = temp;
            if (big == 0.0)
                throw new MatrixException("Matrix is singular");
            vv[i-1] = 1.0/big;
        }
        
        for (j = 1; j <= n; j++) {
            for (i = 1; i < j; i++) {
                sum = a.data[i-1][j-1];
                for (k = 1; k < i; k++)
                    sum -= a.data[i-1][k-1]*a.data[k-1][j-1];
                a.data[i-1][j-1] = sum;
            }
            big = 0.0;
            for (i = j; i <= n; i++) {
                sum = a.data[i-1][j-1];
                for (k = 1; k < j; k++)
                    sum -= a.data[i-1][k-1]*a.data[k-1][j-1];
                a.data[i-1][j-1] = sum;
                if ((dum = vv[i-1]*Math.abs(sum)) >= big) {
                    big  = dum;
                    imax = i;
                }
            }
            if (j != imax) {
                for (k = 1; k <= n; k++) {
                    dum = a.data[imax-1][k-1];
                    a.data[imax-1][k-1] = a.data[j-1][k-1];
                    a.data[j-1][k-1] = dum;
                }
                d[0] = -d[0];
                vv[imax-1] = vv[j-1];
            }
            if (a.data[j-1][j-1] == 0.0)
                a.data[j-1][j-1] = 1.0e-20;
            if (j != n) {
                dum = 1.0/a.data[j-1][j-1];
                for (i = j+1; i <= n; i++)
                    a.data[i-1][j-1] *= dum;
            }
        }
        
        return a;
    }

// Tester function
    public static void main(String[] args) {
	Matrix A = new GeneralMatrix(2,2);
	A.setIJ(0,0,4);
	A.setIJ(0,1,3);
	A.setIJ(1,0,2);
	A.setIJ(1,1,1);
	
	Matrix B = new GeneralMatrix(2,2);
	B.setIJ(0,0,5);
	B.setIJ(0,1,-3);
	B.setIJ(1,0,0);
	B.setIJ(1,1,11);
	
	Matrix ACopy = new GeneralMatrix((GeneralMatrix)A);
	System.out.println("Get 0,0 value of matrix A: "+A.getIJ(0,0)+"\n");
	System.out.println("Matrix A: \n"+A.toString());
	System.out.println("Determinant of A is: "+A.determinant()+"\n");

	System.out.println("Matrix B: \n"+B.toString());
	System.out.println("Determinant of B is: "+B.determinant()+"\n");

	System.out.println("Matrix ACopy: \n"+ACopy.toString());
	
	Matrix C = A.add(B);
	System.out.println("Matrix A+B: \n"+C.toString());

	Matrix D = A.multiply(5);
	System.out.println("Matrix A*5: \n"+D.toString());

	Matrix E = A.multiply(B);
	System.out.println("Matrix A*B: \n"+E.toString());
  
	Matrix F = new GeneralMatrix(3,3);
	F.random();
	System.out.println("Random Matrix: \n"+F.toString());

	System.out.println("Trying to premultiply 2x2 matrix with 3x10 matrix: ");
	try{
	    Matrix a = new GeneralMatrix(3,10);
	    a.random();
	    Matrix b = new GeneralMatrix(2,2);
	    b.random();
	    Matrix c = b.multiply(a);
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}

	System.out.println("\nTrying to add 5x5 matrix to 3x3 matrix: ");
	try{
	    Matrix a1 = new GeneralMatrix(5,5);
	    a1.random();
	    Matrix b1 = new GeneralMatrix(3,3);
	    b1.random();
	    Matrix c1 = b1.add(a1);
	    }catch(MatrixException except){
		System.out.println(except.getMessage());
	    }		

	System.out.println("\nTrying to find determinant of 7x4 matrix: ");
	try{
	    Matrix a2 = new GeneralMatrix(7,4);
	    a2.determinant();
	    }catch(MatrixException except){
		System.out.println(except.getMessage());
	    }	

	System.out.println("\nTrying to find determinant of singular 2x2 matrix: ");
	try{
	    Matrix a3 = new GeneralMatrix(2,2);
	    a3.setIJ(0,0,1);
	    a3.setIJ(0,1,1);
	    a3.setIJ(1,0,0);
	    a3.setIJ(1,1,0);		
	    a3.determinant();
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}	

	System.out.println("\nTrying to Create Matrix with one dimension zero ");
	try{
	    Matrix a4 = new GeneralMatrix(0,7);	
	}catch(MatrixException except){
	    System.out.println(except.getMessage());
	}		  
  }
}
