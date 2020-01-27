
# Subgroup Discovery
# Copyright (C) 2020  Jurian Baas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Credit scoring data.
#'
#' A dataset containing the attributes of people who did or did not default on their loan.
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
"credit"

#' Pima Indians Diabetes Database.
#'
#' Sources:
#' (a) Original owners: National Institute of Diabetes and Digestive and Kidney Diseases
#' (b) Donor of database: Vincent Sigillito \email{vgs@aplcen.apl.jhu.edu}
#'     Research Center, RMI Group Leader
#'     Applied Physics Laboratory
#'     The Johns Hopkins University
#'     Johns Hopkins Road
#'     Laurel, MD 20707
#'     (301) 953-6231
#' (c) Date received: 9 May 1990
#'
#' Past Usage:
#'  Smith, J.W., Everhart, J.E., Dickson, W.C., Knowler, W.C., &
#'  Johannes, R.S. (1988). Using the ADAP learning algorithm to forecast
#'  the onset of diabetes mellitus.  In Proceedings of the Symposium
#'  on Computer Applications and Medical Care (pp. 261--265).
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



#' Ames Housing data.
#'
#' Data set contains information from the Ames Assessor Office used in computing assessed
#' values for individual residential properties sold in Ames, IA from 2006 to 2010.
#' Tab characters are used to separate variables in the data file. The data has 82 columns
#' which include 23 nominal, 23 ordinal, 14 discrete, and 20 continuous variables
#' (and 2 additional observation identifiers).
#'
#' Sources:
#' Ames, Iowa Assessor Office
#'
#' @format A data frame with 2930 rows and 82 variables:
#' \describe{
#'   \item{Order}{(Discrete): Observation number}
#'   \item{PID}{(Nominal): Parcel identification number  - can be used with city web site for parcel review}
#'   \item{MS.SubClass}{(Nominal): Identifies the type of dwelling involved in the sale}
#'   \item{MS.Zoning}{(Nominal): Identifies the general zoning classification of the sale}
#'   \item{Lot.Frontage}{(Continuous): Linear feet of street connected to property}
#'   \item{Lot.Area}{(Continuous): Lot size in square feet}
#'   \item{Street}{(Nominal): Type of road access to property}
#'   \item{Alley}{(Nominal): Type of alley access to property}
#'   \item{Lot.Shape}{(Ordinal): General shape of property}
#'   \item{Land.Contour}{(Nominal): Flatness of the property}
#'   \item{Utilities}{(Ordinal): Type of utilities available}
#'   \item{Lot.Config}{(Nominal): Lot configuration}
#'   \item{Land.Slope}{(Ordinal): Slope of property}
#'   \item{Neighborhood}{(Nominal): Physical locations within Ames city limits (map available)}
#'   \item{Condition.1}{(Nominal): Proximity to various conditions}
#'   \item{Condition.2}{(Nominal): Proximity to various conditions (if more than one is present)}
#'   \item{Bldg.Type}{(Nominal): Type of dwelling}
#'   \item{House.Style}{(Nominal): Style of dwelling}
#'   \item{Overall.Qual}{(Ordinal): Rates the overall material and finish of the house}
#'   \item{Overall.Cond}{(Ordinal): Rates the overall condition of the house}
#'   \item{Year.Built}{(Discrete): Original construction date}
#'   \item{Year.Remod.Add}{(Discrete): Remodel date (same as construction date if no remodeling or additions)}
#'   \item{Roof.Style}{(Nominal): Type of roof}
#'   \item{Roof.Matl}{(Nominal): Roof material}
#'   \item{Exterior.1st}{(Nominal): Exterior covering on house}
#'   \item{Exterior.2nd}{(Nominal): Exterior covering on house (if more than one material)}
#'   \item{Mas.Vnr.Type}{(Nominal): Masonry veneer type}
#'   \item{Mas.Vnr.Area}{(Continuous): Masonry veneer area in square feet}
#'   \item{Exter.Qual}{(Ordinal): Evaluates the quality of the material on the exterior}
#'   \item{Exter.Cond}{(Ordinal): Evaluates the present condition of the material on the exterior}
#'   \item{Foundation}{(Nominal): Type of foundation}
#'   \item{Bsmt.Qual}{(Ordinal): Evaluates the height of the basement}
#'   \item{Bsmt.Cond}{(Ordinal): Evaluates the general condition of the basement}
#'   \item{Bsmt.Exposure}{	(Ordinal): Refers to walkout or garden level walls}
#'   \item{BsmtFin.Type.1}{	(Ordinal): Rating of basement finished area}
#'   \item{BsmtFin.SF.1}{(Continuous): Type 1 finished square feet}
#'   \item{BsmtFin.Type.2}{	(Ordinal): Rating of basement finished area (if multiple types)}
#'   \item{BsmtFin.SF.2}{(Continuous): Type 2 finished square feet}
#'   \item{Bsmt.Unf.SF}{(Continuous): Unfinished square feet of basement area}
#'   \item{Total.Bsmt.SF}{(Continuous): Total square feet of basement area}
#'   \item{Heating}{	(Nominal): Type of heating}
#'   \item{Heating.QC}{(Ordinal): Heating quality and condition}
#'   \item{Central.Air}{(Nominal): Central air conditioning}
#'   \item{Electrical}{(Ordinal): Electrical system}
#'   \item{X1st.Flr.SF}{(Continuous): First Floor square feet}
#'   \item{X2nd.Flr.SF}{(Continuous)	: Second floor square feet}
#'   \item{Low.Qual.Fin.SF}{(Continuous): Low quality finished square feet (all floors)}
#'   \item{Gr.Liv.Area}{(Continuous): Above grade (ground) living area square feet}
#'   \item{Bsmt.Full.Bath}{(Discrete): Basement full bathrooms}
#'   \item{Bsmt.Half.Bath}{(Discrete): Basement half bathrooms}
#'   \item{Full.Bath}{(Discrete): Full bathrooms above grade}
#'   \item{Half.Bath}{(Discrete): Half baths above grade}
#'   \item{Bedroom.AbvGr}{(Discrete): Bedrooms above grade (does NOT include basement bedrooms)}
#'   \item{Kitchen.AbvGr}{(Discrete): Kitchens above grade}
#'   \item{Kitchen.Qual}{(Ordinal): Kitchen quality}
#'   \item{TotRms.AbvGrd}{	(Discrete): Total rooms above grade (does not include bathrooms)}
#'   \item{Functional}{(Ordinal): Home functionality (Assume typical unless deductions are warranted)}
#'   \item{Fireplaces}{(Discrete): Number of fireplaces}
#'   \item{Fireplace.Qu}{(Ordinal): Fireplace quality}
#'   \item{Garage.Type}{(Nominal): Garage location}
#'   \item{Garage.Yr.Blt}{(Discrete): Year garage was built}
#'   \item{Garage.Finish}{(Ordinal)	: Interior finish of the garage}
#'   \item{Garage.Cars}{(Discrete): Size of garage in car capacity}
#'   \item{Garage.Area}{(Continuous): Size of garage in square feet}
#'   \item{Garage.Qual}{(Ordinal): Garage quality}
#'   \item{Garage.Cond}{(Ordinal): Garage condition}
#'   \item{Paved.Drive}{(Ordinal): Paved driveway}
#'   \item{Wood.Deck.SF}{(Continuous): Wood deck area in square feet}
#'   \item{Open.Porch.SF}{(Continuous): Open porch area in square feet}
#'   \item{Enclosed.Porch}{(Continuous): Enclosed porch area in square feet}
#'   \item{X3Ssn.Porch}{(Continuous): Three season porch area in square feet}
#'   \item{Screen.Porch}{(Continuous): Screen porch area in square feet}
#'   \item{Pool.Area}{(Continuous): Pool area in square feet}
#'   \item{Pool.QC}{(Ordinal): Pool quality}
#'   \item{Fence}{(Ordinal): Fence quality}
#'   \item{Misc.Feature}{(Nominal): Miscellaneous feature not covered in other categories}
#'   \item{Misc.Val}{(Continuous): $Value of miscellaneous feature}
#'   \item{Mo.Sold}{(Discrete): Month Sold (MM)}
#'   \item{Yr.Sold}{(Discrete): Year Sold (YYYY)}
#'   \item{Sale.Type}{(Nominal): Type of sale}
#'   \item{Sale.Condition}{(Nominal): Condition of sale}
#'   \item{SalePrice}{(Continuous): Sale price}
#' }
"ames"
