\name{CoSIAn-class}
\Rdversion{1.1}
\docType{class}
\alias{CoSIAn-class}
\alias{getConversion,CoSIAn-method}
\alias{getGEx,CoSIAn-method}
\alias{getGExMetrics,CoSIAn-method}
\alias{plotCVGEx,CoSIAn-method}
\alias{plotDSGEx,CoSIAn-method}
\alias{plotSpeciesGEx,CoSIAn-method}
\alias{plotTissueGEx,CoSIAn-method}

\title{Class \code{"CoSIAn"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CoSIAn", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{gene_set}:}{Object of class \code{"character"} ~~ }
    \item{\code{i_species}:}{Object of class \code{"character"} ~~ }
    \item{\code{input_id}:}{Object of class \code{"character"} ~~ }
    \item{\code{o_species}:}{Object of class \code{"character"} ~~ }
    \item{\code{output_ids}:}{Object of class \code{"character"} ~~ }
    \item{\code{mapping_tool}:}{Object of class \code{"character"} ~~ }
    \item{\code{ortholog_database}:}{Object of class \code{"character"} ~~ }
    \item{\code{converted_id}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{map_tissues}:}{Object of class \code{"character"} ~~ }
    \item{\code{map_species}:}{Object of class \code{"character"} ~~ }
    \item{\code{gex}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{metric_type}:}{Object of class \code{"character"} ~~ }
    \item{\code{metric}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CoSIA}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getConversion}{\code{signature(object = "CoSIAn")}: ... }
    \item{getGEx}{\code{signature(object = "CoSIAn")}: ... }
    \item{getGExMetrics}{\code{signature(object = "CoSIAn")}: ... }
    \item{plotCVGEx}{\code{signature(object = "CoSIAn")}: ... }
    \item{plotDSGEx}{\code{signature(object = "CoSIAn")}: ... }
    \item{plotSpeciesGEx}{\code{signature(object = "CoSIAn")}: ... }
    \item{plotTissueGEx}{\code{signature(object = "CoSIAn")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("CoSIAn")
}
\keyword{classes}
