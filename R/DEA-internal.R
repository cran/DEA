`.check` <-
function (X , Y) {

	if ( length(dim(X)) != 2 ) stop ("Non-bidimensional inputs matrix")

	if ( length(dim(Y)) != 2 ) stop ("Non-bidimensional outputs matrix")

	if ( dim(X)[1] != dim(Y)[1] ) stop ("The matrices of inputs and outputs have a different number of rows")

	if (any(row.names(X) != row.names(Y))) warning ("The matrices of inputs and outputs have a different name of DMUs")

	if ( any ( is.na(X) ) ) stop ("Inputs matrix contains NAs")


	if ( any ( is.na(Y) ) ) stop ("Outputs matrix contains NAs")
}

`.primera` <-
function (X) {
  if (! is.data.frame(X)) {
    X <- data.frame(X);
    names(X) <- paste("x",1:dim(X)[2], sep="");
    row.names(X) <- paste ("DMU", 1:dim(X)[1] , sep ="" );
  }
  X
}

`.segunda` <-
function (Y) {
  if (! is.data.frame(Y)) {
    Y <- data.frame(Y);
    names(Y) <- paste("y",1:dim(Y)[2], sep="");
    row.names(Y) <- paste ("DMU", 1:dim(Y)[1] , sep ="" );
  }
  Y
}

