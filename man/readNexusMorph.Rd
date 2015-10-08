\name{readNexusMorph}
\alias{readNexusMorph}
\title{Reads morphological character data In NEXUS Format}
\description{This function reads a file with character scoring in the NEXUS format and provides a properly formatted data set to input into the lewisMkv function.}
\usage{readNexusMorph(file)}
\arguments{
\item{file}{a file name specified by either a variable of mode character, or a double-quoted string}
}
\details{
This parser is really just a hacked version of \dQuote{read.nexus.data} that is provided by \code{ape}. I rewrote parts to deal with completely and incompletely ambiguous characters, gaps, and so that the format of the data fits perfectly into \code{lewisMkv}. As with \dQuote{read.nexus.data} this parser reads data from a file written in a \emph{restricted} NEXUS format. I've posted on my website (http://www.jeremybeaulieu.org/r.html) an example NEXUS file properly formatted that will work with this function. But there is one major change worth mentioning here:
\itemize{
\item{Multistate characters ARE allowed. That is, NEXUS allows you to specify multiple character states at a character position either as an uncertainty or as an actual appearance of multiple states, \code{\{XY\}}. Examples:\cr
\code{taxon 0011?110} \bold{--- OK!}\cr
\code{taxon 0011-110} \bold{--- OK!}\cr
\code{taxon 0011{01}110} \bold{--- OK!}\cr
\code{taxon 00?1{01}1-0} \bold{--- OK!}\cr
\code{taxon 0011(01)110} \bold{--- STILL NOT OK!}}
}
}
\value{
A data.frame of character scorings across sites for a set of species. Gaps and completely missing characters are coded as \dQuote{?}, and \dQuote{&} separate scorings of partially ambiguous characters. Will work as a direct input for \dQuote{lewisMkv}.
}
\author{Jeremy M. Beaulieu}
\keyword{file}
