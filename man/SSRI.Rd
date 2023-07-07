\name{SSRI}
\alias{SSRI}
\docType{data}
\title{
Antidepressants and bone density
}
\description{
NHANES 2009-2010 data on selective serotonin reuptake inhibitors (SSRIs) and bone density for adults.
}
\usage{data("SSRI")}
\format{
  A data frame with 2759 observations on the following 11 variables.
  \describe{
    \item{\code{z}}{Antidepressants status, 1 = SSRI user, 0 = taking neither SSRI nor TCA}
    \item{\code{female}}{1 = female, 0 = male}
    \item{\code{age}}{Age in years, >=18}
    \item{\code{black}}{1=black race, 0=other}
    \item{\code{hispanic}}{1=Hispanic race, 0=other}
    \item{\code{povertyNA}}{1=missing poverty ratio, 0=observed poverty ratio}
    \item{\code{povertyFill}}{Ratio of family income to the poverty level, capped at 5 times poverty, with missing values imputed as mean}
    \item{\code{education}}{Level of education in 5-point scale, with 1 for less than 9th grade, 3 for high school or equivalent, and 5 for college graduate}
    \item{\code{height}}{Height in cm}
    \item{\code{weight}}{Weight in kg}
    \item{\code{bmi}}{BMI or body-mass-index}
    \item{\code{cotinine}}{serum cotinine level (ng/mL)}
    \item{\code{hd}}{serum 25-hydroxyvitamin D2 + D3 level (nmol/L)}
    \item{\code{diabetes}}{1=diagnosed with diabetes or borderline diabetes, 0=no diabetes}
    \item{\code{insurance}}{whether covered by insurance}
    \item{\code{weighmore}}{whether would like to weigh more}
    \item{\code{weighless}}{whether would like to weigh less}
    \item{\code{weightchange}}{weight change in the recent 1 year of finishing the questionnaire (kg)}
    \item{\code{physicalact}}{whether has moderate and vigorous activity in the past 30 days}
    \item{\code{dietsup}}{whether has taken any dietary supplements in the past month}
    \item{\code{pr}}{Height}
    \item{\code{femurbmd}}{Femur bone mineral density}
    \item{\code{femurbmc}}{Femur bone mineral content}
    \item{\code{fneckbmd}}{Femoral neck bone mineral density}
    \item{\code{fneckbmc}}{Femoral neck bone mineral content}
  }
}

\source{
From the NHANES web page, for NHANES 2009-2010.
}
\references{
US National Health and Nutrition Examination Survey, 2009-2010.  From the US Center for Health Statistics.
}
\examples{
data(SSRI)
summary(SSRI)
}
\keyword{datasets}
