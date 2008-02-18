`dea.ccr.io.env` <-
function (X, Y, pslv = FALSE , dual = FALSE , infor = FALSE ) {

  X <- .primera (X);
  Y <- .segunda (Y);
  .check(X,Y);
  
  pslv <- as.logical (pslv );
  presolver <- 0;
  if ( pslv ) presolver <- 1;

  dual <- as.logical ( dual );
  dual.simplex <- 0;
  if ( dual ) dual.simplex <- 1;
  
  infor <- as.logical ( infor );
  information <- 0;
  if ( infor ) information <- 1;
  
  nombre.DMUs <- row.names(X);
  nombre.inputs <- names(X);
  nombre.outputs <- names(Y);

  num.DMUs <- dim(X)[1];
  num.inputs <- dim(X)[2];
  num.outputs <- dim(Y)[2];

  XX <- as.matrix(X);
  YY <- as.matrix(Y);
  ZZ <- double(num.DMUs * (num.DMUs + 2));
  

  resul <- .C( "ccr_io_env",
              as.character (nombre.DMUs),
              as.integer (num.DMUs),
              as.character (nombre.inputs),
              as.integer (num.inputs),
              as.character(nombre.outputs),
              as.integer (num.outputs),
              as.double (XX),
              as.double (YY),
	      as.integer (presolver),
              as.integer (dual.simplex),
              as.integer (information),
              salida = as.double (ZZ),
	      PACKAGE = "DEA");

  resul <- resul$salida;
  resul <- matrix(resul, num.DMUs , num.DMUs + 2 );

  lambda <- data.frame( resul[ , 1:num.DMUs ] );
  eff <- resul[ , num.DMUs+1 ];
  slack <- resul[ , num.DMUs+2 ];
  
  row.names(lambda) <- nombre.DMUs;
  nLs <- paste ( rep ("L.", num.DMUs) , nombre.DMUs , sep = "");
  names(lambda) <- nLs;

  names(eff) <- nombre.DMUs;
  names(slack) <- nombre.DMUs;
  list ( eff = eff , lambda = lambda , slack = slack );
  
}

