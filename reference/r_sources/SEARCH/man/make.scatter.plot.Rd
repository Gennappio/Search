\name{make.scatter.plot}
\alias{make.scatter.plot}
\title{
make.scatter.plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.scatter.plot(x = NULL, y = NULL, correlation.method = "spearman", correlation.coef = TRUE, regression.line = TRUE, lowess.line = TRUE, xlimit = NULL, x.abline.at = NULL, y.abline.at = NULL, cex.axis = 1.2, cex.lab = 1.4, main.title = "", xaxis.label = expression(bold(paste('log'[2], ' mRNA abundance'))), yaxis.label = expression(bold(paste('log'[2], ' mRNA abundance'))), yaxis.label.offset = -1.4, margins = c(2.7, 3.5, 0.5, 0.2), file.name = NULL, height = 5, width = 5, resolution = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of fold change in log space, ideally log2}
  \item{y}{A vector of P values or adjusted P values in -log10 space}
  \item{correlation.method}{A string specifying which correlation method to use. Supported methods are same as in \code{cor.test}. Defaults to spearman}
  \item{correlation.coef}{Logical indicating whether to display correlation coefficient on the plot}
  \item{regression.line}{Logical indicating whether to fit and plot regression line to the data}
  \item{lowess.line}{Logical indicating whether to fit and plot Lowess line to the data}
  \item{xlimit}{A vector of length two indicate the min and max limit of x-axis}
  \item{x.abline.at}{
%%     ~~Describe \code{x.abline.at} here~~
}
  \item{y.abline.at}{
%%     ~~Describe \code{y.abline.at} here~~
}
  \item{cex.axis}{
%%     ~~Describe \code{cex.axis} here~~
}
  \item{cex.lab}{
%%     ~~Describe \code{cex.lab} here~~
}
  \item{main.title}{
%%     ~~Describe \code{main.title} here~~
}
  \item{xaxis.label}{
%%     ~~Describe \code{xaxis.label} here~~
}
  \item{yaxis.label}{
%%     ~~Describe \code{yaxis.label} here~~
}
  \item{yaxis.label.offset}{}
  \item{margins}{A vector of length four indicating space margin: bottom, left, top, right}
  \item{file.name}{A string specifying the file name}
  \item{height}{A number specifying height of the plot in inches}
  \item{width}{A number specifying width of the plot in inches}
  \item{resolution}{A number specifying resolution of the plot}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Syed Haider
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
