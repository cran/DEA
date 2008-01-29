\name{dea.add.mul}
\alias{dea.add.mul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Additive model in multiplier form }
\description{
  Fits an additive DEA model in multiplier form. This model corresponds to eqs. 4.39-4.43 in Cooper \emph{et al}., 2007 (see References).
}
\usage{
dea.add.mul(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ inputs. Bidimensional data frame or object that can be coerced to a data frame. Each row corresponds to a DMU, each column to an input. }
  \item{Y}{ outputs. Bidimensional data frame or object that can be coerced to a data frame. Each row corresponds to a DMU, each column to an output. }
}
\details{
  Missing values are not allowed.
}
\value{
  a list with components:
  \item{weights }{a data frame containing the weights. Each row corresponds to a DMU and each column to a weight.}
  \item{slack}{a vector containing the sum of slacks for each DMU.}
}
\references{

	Cooper, W. W., Seiford, L. M. and Tone, K.  (2007):
	\emph{Data envelopment analysis: a comprehensive text with models, applications, references and DEA-solver software.} Second edition. Springer. New York.
}
\seealso{\code{\link{dea.add.env}}}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
