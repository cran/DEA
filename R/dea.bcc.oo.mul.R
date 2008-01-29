`dea.bcc.oo.mul` <-
function (X, Y) {
  
  X <- .primera (X);
  Y <- .segunda (Y);
  .check(X,Y);
  
  nombre.DMUs <- row.names(X);
  nombre.inputs <- names(X);
  nombre.outputs <- names(Y);

  num.DMUs <- dim(X)[1];
  num.inputs <- dim(X)[2];
  num.outputs <- dim(Y)[2];

  ll <- num.inputs + num.outputs + 1  

  XX <- as.matrix(X);
  YY <- as.matrix(Y);
  ZZ <- double(num.DMUs * ll );

  resul <- .C( "bcc_oo_mul",
              as.character (nombre.DMUs),
              as.integer (num.DMUs),
              as.character (nombre.inputs),
              as.integer (num.inputs),
              as.character(nombre.outputs),
              as.integer (num.outputs),
              as.double (XX),
              as.double (YY),
              salida = as.double (ZZ),
	      PACKAGE = "DEA");

  resul <- resul$salida;
  resul <- matrix(resul, num.DMUs , ll );

  weights <- data.frame( resul[ , -ll ] )
  eff <- resul[, ll ]
  
  
  row.names(weights) <- nombre.DMUs;
  names( eff ) <- nombre.DMUs;

  nin  <- paste (rep("v.", num.inputs ) , nombre.inputs , sep ="" );
  nout <- paste (rep("u.", num.outputs) , nombre.outputs, sep ="" );

  names(weights) <- c( nin , nout );
  list( eff = eff , weights = weights );
  
}

