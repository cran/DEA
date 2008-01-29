.First.lib <- function ( ll , pp ) {

	library.dynam("DEA" , pp , ll)

}

.Last.lib <- function ( ruta ) {

	library.dynam.unload("DEA" , ruta )

}
