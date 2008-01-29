\name{dea.bcc.oo.mul}
\alias{dea.bcc.oo.mul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BCC output oriented model in multiplier form }
\description{
  Fits an output oriented DEA model in multiplier form with variable returns to scale. This model corresponds to eqs. 4.27-4.30 in Cooper \emph{et al}., 2007 (see References).
}
\usage{
dea.bcc.oo.mul(X, Y)
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
  \item{eff }{a vector containing the technical efficiency of the DMUs.}
  \item{weights }{a data frame containing the weights. Each row corresponds to a DMU and each column to a weight.}
}
\references{

	Cooper, W. W., Seiford, L. M. and Tone, K.  (2007):
	\emph{Data envelopment analysis: a comprehensive text with models, applications, references and DEA-solver software.} Second edition. Springer. New York.
}

\seealso{ \code{\link{dea.bcc.oo.env}}}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
