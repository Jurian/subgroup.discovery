#' Credit scoring data.
#'
#' A dataset containing the atributes of people who did or did not default on their loan.
#'
#' Toy example, useful for debugging purposes.
#'
#' @format A data frame with 10 rows and 6 variables:
#' \describe{
#'   \item{age}{age of person}
#'   \item{married}{is the person married or not, boolean}
#'   \item{house}{is the person a homeowner or not, boolean}
#'   \item{income}{income of person, in thousands}
#'   \item{gender}{factor with levels male, female}
#'   \item{class}{class variable (0 or 1)}
#' }
#' @source \url{http://www.cs.uu.nl/docs/vakken/adm/bump.pdf}
"credit"

#' Pima Indians Diabetes Database.
#'
#' Sources:
#' (a) Original owners: National Institute of Diabetes and Digestive and Kidney Diseases
#' (b) Donor of database: Vincent Sigillito (vgs@aplcen.apl.jhu.edu)
#'     Research Center, RMI Group Leader
#'     Applied Physics Laboratory
#'     The Johns Hopkins University
#'     Johns Hopkins Road
#'     Laurel, MD 20707
#'     (301) 953-6231
#' (c) Date received: 9 May 1990
#'
#' Past Usage:
#'  Smith,~J.~W., Everhart,~J.~E., Dickson,~W.~C., Knowler,~W.~C., \&
#'  Johannes,~R.~S. (1988). Using the ADAP learning algorithm to forecast
#'  the onset of diabetes mellitus.  In {Proceedings of the Symposium
#'  on Computer Applications and Medical Care} (pp. 261--265).
#'  IEEE Computer Society Press.
#'
#'  The diagnostic, binary-valued variable investigated is whether the
#'  patient shows signs of diabetes according to World Health Organization
#'  criteria (i.e., if the 2 hour post-load plasma glucose was at least
#'  200 mg/dl at any survey  examination or if found during routine medical
#'  care). The population lives near Phoenix, Arizona, USA.
#'
#'  Results: Their ADAP algorithm makes a real-valued prediction between
#'  0 and 1.  This was transformed into a binary decision using a cutoff of
#'  0.448.  Using 576 training instances, the sensitivity and specificity
#'  of their algorithm was 76% on the remaining 192 instances.
#'
#' Relevant Information:
#'  Several constraints were placed on the selection of these instances from
#'  a larger database.  In particular, all patients here are females at
#'  least 21 years old of Pima Indian heritage.
#'
#' @format A data frame with 768 rows and 9 variables:
#' \describe{
#'   \item{pregnant}{Number of times pregnant}
#'   \item{glucose}{Plasma glucose concentration a 2 hours in an oral glucose tolerance test}
#'   \item{bp}{Diastolic blood pressure (mm Hg)}
#'   \item{skin_thickness}{Triceps skin fold thickness (mm)}
#'   \item{insulin}{2-Hour serum insulin (mu U/ml)}
#'   \item{bmi}{Body mass index (weight in kg/(height in m)^2)}
#'   \item{diabetes}{Diabetes pedigree function}
#'   \item{age}{Age (years)}
#'   \item{class}{Class variable (0 or 1)}
#' }
#' @source \url{https://archive.ics.uci.edu/ml/datasets/pima+indians+diabetes}
"pima"
