\name{dea.bcc.oo.env}
\alias{dea.bcc.oo.env}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ BCC output oriented model in envelopment form }
\description{
  Fits an output oriented DEA model in envelopment form with variable returns to scale. This model corresponds to eqs. 4.22-4.26 in Cooper \emph{et al}., 2007 (see References).
}
\usage{
dea.bcc.oo.env( X , Y , pslv = FALSE , dual = FALSE , infor = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ inputs. Bidimensional data frame or object that can be coerced to a data frame. Each row corresponds to a DMU, each column to an input. }
  \item{Y}{ outputs. Bidimensional data frame or object that can be coerced to a data frame. Each row corresponds to a DMU, each column to an output. }
  \item{pslv}{ logical. Defaults to \code{FALSE}. If \code{TRUE}, the simplex presolver will be used.}
  \item{dual}{ logical. Defaults to \code{FALSE}. IF \code{TRUE} and the initial basic solution is dual feasible, dual simplex is used. }
  \item{infor}{ logical. Defaults to \code{FALSE}. If \code{TRUE}, information about the convergence of the optimization process will be sent to the console. }
}
\details{
  Missing values are not allowed.\newline \newline If \code{pslv=TRUE}, a glpk built-in LP presolver is enabled. It may improve the numeric properties of the LP problem. \newline \newline
  If \code{infor=TRUE}, each line of information sent to the console consists of ten columns with some details about the calculus of the efficiency for each DMU. First five columns correspond to the phase I of the linear programming problem and the rest correspond to the phase II (see Cooper \emph{et al}., 2007, sec. 3.3). First and sixth columns indicate if the linear programs have been successfully solved, or not ("OK" or "FAULT"). The other ones contain information reported by the glpk solver about the quality of the solutions. Their values can be "H", "M", "L" or "?", which mean high, medium, low or wrong, respectively. If all the indicators show high or medium quality, the user can be sure that the solution is quite accurate. For more details see Makhorin (2005).
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

	Makhorin, A. (2005): \emph{GNU Linear Programming Kit. Reference Manual. Version 4.8.} Available at \url{http://www.gnu.org/software/glpk/}
}

\seealso{ \code{\link{dea.bcc.oo.mul}}}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
