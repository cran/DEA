\name{dea.sbm.ccr.io}
\alias{dea.sbm.ccr.io}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ CCR input oriented slack-based measure of efficiency }
\description{
  Fits an input oriented slack-based DEA model with constant returns to scale. This model corresponds to eq. 4.78 in Cooper \emph{et al}., 2007 (see References).
}
\usage{
dea.sbm.ccr.io(X, Y)
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
  \item{ro}{a vector containing the slack-based-measure of efficiency of the DMUs.}
  \item{lambda}{a bidimensional data frame containing the  \eqn{\lambda} coefficients. Each row corresponds to a DMU. The columns correspond to the \eqn{\lambda}s.}
  \item{slack}{a bidimensional data frame containing the slacks for each DMU. Each row corresponds to a DMU.}
}
\references{

	Cooper, W. W., Seiford, L. M. and Tone, K.  (2007):
	\emph{Data envelopment analysis: a comprehensive text with models, applications, references and DEA-solver software.} Second edition. Springer. New York.
}
\seealso{ \code{\link{dea.sbm.bcc.io}}}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
