\name{readNexusMorph}
\alias{readNexusMorph}
\title{
  Reads morphological character data In NEXUS Format
}
\description{
  This function reads a file with character scoring in the NEXUS format and provides a properly formatted data set to input into the lewisMkv function. 
}
\usage{
readNexusMorph(file)
}
\arguments{
\item{file}{a file name specified by either a variable of mode character, or a double-quoted string}
}
\details{
This parser is really just a hacked version of \dQuote{read.nexus.data} that is provided by \code{ape}. I rewrote parts to deal with missing, complete and incomplete ambiguous characters, and so that the format of the data fits perfectly into \code{lewisMkv}. As with \dQuote{read.nexus.data} this parser read data from a file written in a \emph{restricted} NEXUS format. I've included an example NEXUS formatted file that will work with this function on my website: http://www.jeremybeaulieu.org/r.html. But there is one major change worth mentioning here:

\itemize{
\item{\bold{I}}{Multistate characters ARE allowed. That is, NEXUS allows you to specify multiple character states at a character position either as an uncertainty, \code{(XY)}, or as an actual appearance of multiple states, \code{\{XY\}}. Examples:\cr
\code{taxon 0011?110} \bold{--- OK!}\cr
\code{taxon 0011-110} \bold{--- OK!}\cr
\code{taxon 0011{01}110} \bold{--- OK!}\cr
\code{taxon 00?1{01}1-0} \bold{--- OK!}\cr
\code{taxon 0011(01)110} \bold{--- STILL NOT OK!}}
}
}
\value{
A data.frame of character scoring across sites for each species. Gaps and completely missing characters are coded as \dQuote{?}, and partially ambiguous character scorings are separated by \dQuote{&}. Will work as a direct input for \dQuote{lewisMkv}.
}
\author{Jeremy M. Beaulieu}
\keyword{file}
