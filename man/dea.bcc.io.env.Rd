\name{dea.bcc.io.env}
\alias{dea.bcc.io.env}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BCC input oriented model in envelopment form }
\description{
  Fits an input oriented DEA model in envelopment form with variable returns to scale. This model corresponds to eqs. 4.2-4.6 in Cooper \emph{et al}., 2007 (see References).
}
\usage{
dea.bcc.io.env(X, Y)
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
  \item{eff}{a vector containing the technical efficiency of the DMUs.}
  \item{lambda}{a bidimensional data frame containing the  \eqn{\lambda} coefficients. Each row corresponds to a DMU. The columns correspond to the \eqn{\lambda}s.}
  \item{slack}{a vector containing the sum of slacks for each DMU.}
}
\references{

	Cooper, W. W., Seiford, L. M. and Tone, K.  (2007):
	\emph{Data envelopment analysis: a comprehensive text with models, applications, references and DEA-solver software.} Second edition. Springer. New York.
}
\seealso{ \code{\link{dea.bcc.io.mul}}}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
