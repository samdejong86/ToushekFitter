/*
 *  This header contains methods to perform a least squares fit.
 *
 *
 */


//the TMatrix multiplier doesn't seem to work...
TMatrixD MatrixMultiply(TMatrixD a, TMatrixD b){
  
  int r1=a.GetNrows();
  int r2=b.GetNrows();
  int c1=a.GetNcols();
  int c2=b.GetNcols();
  
  TMatrixD mult(r1,c2);
  
  int i,j,k;
  
  /* Initializing elements of matrix mult to 0.*/
  for(i=0; i<r1; ++i){
    for(j=0; j<c2; ++j){
      mult[i][j]=0;
    }
  }
  /* Multiplying matrix a and b and storing in array mult. */
  for(i=0; i<r1; ++i){
    for(j=0; j<c2; ++j){
      for(k=0; k<c1; ++k){
	mult[i][j]+=a[i][k]*b[k][j];
      }
    }
  }
    
  return mult;
}

//apply the least squares fit
TMatrixD fitter(TMatrixD ydata, TMatrixD X){

  TMatrixD Xt(X.GetNcols(), X.GetNrows());
  for(int i=0; i<X.GetNrows(); i++){
    for(int j=0; j<X.GetNcols(); j++){
      Xt[j][i] = X[i][j];
    }
  }

  //multiply X-transpose by X, then invert
  TMatrixD XtX = MatrixMultiply(Xt, X);

  TMatrixD XtXInv = XtX.Invert();

  //matrix multiply X-transpose by ydata
  TMatrixD XtY = MatrixMultiply(Xt, ydata);

  //matrix product of XtX and XtY is the least squares solution to the fit
  TMatrixD soln = MatrixMultiply(XtXInv, XtY);

  return soln;
  
}



//get the uncertainty
TMatrixD uncertainty(TMatrixD ydata, TMatrixD appliedSoln, TMatrixD X){
  
  TMatrixD Xt(X.GetNcols(), X.GetNrows());
  for(int i=0; i<X.GetNrows(); i++){
    for(int j=0; j<X.GetNcols(); j++){
      Xt[j][i] = X[i][j];
    }
  }

  //multiply X-transpose by X, then invert
  TMatrixD XtX = MatrixMultiply(Xt, X);
  TMatrixD XtXInv = XtX.Invert();

  //get the fitting uncertainties  
  TMatrixD difference = ydata-appliedSoln;

  double sum;
  double sumSq;
  int n=difference.GetNrows();
  for(int i=0; i<n; i++){
    sum +=difference[i][0];
    sumSq+=difference[i][0]*difference[i][0];      
  }
  
  double sigmaSq = (sumSq/n)-(sum/n)*(sum/n);

  //off diagonal terms are the covariances
  TMatrixD Var = XtXInv*sigmaSq;

  return Var;

}
