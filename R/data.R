#' @docType data
#' @author Yujiao Mai, Deo Kumar Srivastava, and Kevin R. Krull
#' @title PedsQL Multidimensional Fatigue Scale Factor Structure
#' @description
#' The data provide the information needed for estimating the CR coefficient Omega-generic of the PedsQL Multidimensional Fatigue Scale \cite{(Varni et al., 2002)}.
#' The estimated parameter matrices (\code{Lambda}, \code{Phi}, and \code{Psi}) were obtained by fitting factor models with participants' responses to the PedsQL Multidimensional Fatigue Scale.
#' Two different factor structures, a three-correlated-factor model and a bi-factor model, were included in the example.
#' Exploratory structural equation modeling \cite{(ESEM; Asparouhov, & Muthen, 2009; Morin, Arens, & Marsh, 2016)} was employed to estimate the model.
#' The sample included 87 young-adult cancer survivors. Sample data were collected by St. Jude LIFE Study \cite{(SJCRH., 2007-2021)}.
#' Please refer to the publication \cite{(Mai, Srivastava, & Krull, 2021)} for more information.
#'
#'
#' @format \code{PedsQLMFS}: A list including three sub-lists: \code{ScaleStructure}, \code{ESEM}, and \code{biESEM}.
#'
#'  1. \code{PedsQLMFS$ScaleStructure}:
#'    \code{ScaleStructure} is a list used to describe the subscale names and items within each subscale. It contains three vectors: \code{GeneralFatigue}, \code{SleepFatigue} and \code{CognitiveFatigue}.
#'      \describe{
#'      \item{\code{GeneralFatigue}}{A vector of item varibale names that are in the subscale "General Fatigue"}
#'      \item{\code{SleepFatigue}}{A vector of item varibale names that are in the subscale "Sleep/rest Fatigue"}
#'      \item{\code{CognitiveFatigue}}{A vector of item varibale names that are in the subscale "Cognitive Fatigue"}
#'      }
#'
#'
#'  2. \code{PedsQLMFS$ESEM}:
#'     \code{ESEM} is a list of parameter matrices of a three-correlated-factor model. It contains three matrices: \code{Lambda}, \code{Phi}, and \code{Psi}.
#'
#'    \code{Lambda}: The factor-loading matrix; A matrix with 18 rows and 3 columns, each row represent one scale item, each column represent one factor.
#'      \describe{
#'      \item{\code{GeneralFatigue}}{Factor loadings on the sub-domain construct "Gneral Fatigue"}
#'      \item{\code{SleepFatigue}}{Factor loadings on the sub-domain construct "Sleep/rest Fatigue"}
#'      \item{\code{CognitiveFatigue}}{Factor loadings on the sub-domain construct "Cognitive Fatigue"}
#'      }
#'
#'    \code{Phi}: The factor variance-covariance matrix; A matrix with 3 rows and 3 columns. Each row represent one factor. So does each column.
#'      \describe{
#'      \item{\code{GlobalFatigue}}{Factor loadings on the global (general factor) construct "Global Fatigue"}
#'      \item{\code{GeneralFatigue}}{Factor loadings on the specific (group factor) construct "General Fatigue"}
#'      \item{\code{SleepFatigue}}{Factor loadings on the specific (group factor) construct "Sleep/rest Fatigue"}
#'      \item{\code{CognitiveFatigue}}{Factor loadings on the specific (group factor) construct "Cognitive Fatigue"}
#'      }
#'
#'    \code{Psi}: The item-error variane-covariance matrix; A matrix with 18 rows and 18 columns. Each row represent one item. So does each column.
#'      \describe{
#'      \item{\code{Y1}}{item GeneralFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y2}}{item GeneralFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y3}}{item GeneralFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y4}}{item GeneralFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y5}}{item GeneralFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y6}}{item GeneralFatigue6 measurement-error variance and covariance with other items}
#'      \item{\code{Y7}}{item SleepFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y8}}{item SleepFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y9}}{item SleepFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y10}}{item SleepFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y11}}{item SleepFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y12}}{item SleepFatigue6 measurement-error variance and covariance with other items}
#'      \item{\code{Y13}}{item CognitiveFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y14}}{item CognitiveFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y15}}{item CognitiveFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y16}}{item CognitiveFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y17}}{item CognitiveFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y18}}{item CognitiveFatigue6 measurement-error variance and covariance with other items}
#'      }
#'
#'
#'  3. \code{PedsQLMFS$biESEM}:
#'     \code{biESEM} is a list of parameter matrices of a bi-factor model. It contains three matrices: \code{Lambda}, \code{Phi}, and \code{Psi}.
#'
#'    \code{Lambda}: The factor-loading matrix; A matrix with 18 rows and 4 columns, each row represent one scale item, each column represent one factor. The first factor is the global factor (also called general factor) of a bi-factor structure .
#'      \describe{
#'      \item{\code{GlobalFatigue}}{Factor loadings on the global (general factor) construct "Global Fatigue"}
#'      \item{\code{GeneralFatigue}}{Factor loadings on the specific (group factor) construct "Gneral Fatigue"}
#'      \item{\code{SleepFatigue}}{Factor loadings on the specific (group factor) construct "Sleep/rest Fatigue"}
#'      \item{\code{CognitiveFatigue}}{Factor loadings on the specific (group factor) construct "Cognitive Fatigue"}
#'      }
#'
#'    \code{Phi}: The factor variance-covariance matrix; A matrix with 4 rows and 4 columns, each row represent one factor, each column represent one factor. The first factor is the global factor (also called general factor) of a bi-factor structure .
#'      \describe{
#'      \item{\code{GlobalFatigue}}{Factor loadings on the global (general factor) construct "Global Fatigue"}
#'      \item{\code{GeneralFatigue}}{Factor loadings on the specific (group factor) construct "General Fatigue"}
#'      \item{\code{SleepFatigue}}{Factor loadings on the specific (group factor) construct "Sleep/rest Fatigue"}
#'      \item{\code{CognitiveFatigue}}{Factor loadings on the specific (group factor) construct "Cognitive Fatigue"}
#'      }
#'
#'    \code{Psi}: The item-error variane-covariance matrix; A matrix with 18 rows and 18 columns. Each row represent one item. So does each column.
#'      \describe{
#'      \item{\code{Y1}}{item GeneralFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y2}}{item GeneralFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y3}}{item GeneralFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y4}}{item GeneralFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y5}}{item GeneralFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y6}}{item GeneralFatigue6 measurement-error variance and covariance with other items}
#'      \item{\code{Y7}}{item SleepFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y8}}{item SleepFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y9}}{item SleepFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y10}}{item SleepFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y11}}{item SleepFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y12}}{item SleepFatigue6 measurement-error variance and covariance with other items}
#'      \item{\code{Y13}}{item CognitiveFatigue1 measurement-error variance and covariance with other items}
#'      \item{\code{Y14}}{item CognitiveFatigue2 measurement-error variance and covariance with other items}
#'      \item{\code{Y15}}{item CognitiveFatigue3 measurement-error variance and covariance with other items}
#'      \item{\code{Y16}}{item CognitiveFatigue4 measurement-error variance and covariance with other items}
#'      \item{\code{Y17}}{item CognitiveFatigue5 measurement-error variance and covariance with other items}
#'      \item{\code{Y18}}{item CognitiveFatigue6 measurement-error variance and covariance with other items}
#'      }
#'
#' @examples
#'
#'  OmegaG::PedsQLMFS$ScaleStructure
#' #  $GeneralFatigue
#' #  [1] "Y1" "Y2" "Y3" "Y4" "Y5" "Y6"
#' #
#' #  $SleepFatigue
#' #  [1] "Y7"  "Y8"  "Y9"  "Y10" "Y11" "Y12"
#' #
#' #  $CognitiveFatigue
#' #  [1] "Y13" "Y14" "Y15" "Y16" "Y17" "Y18"
#'
#'  OmegaG::PedsQLMFS$ESEM$Lambda
#' #           GeneralFatigue  SleepFatigue  CognitiveFatigue
#' #  Y1           0.582        0.134           -0.093
#' #  Y2           0.640        0.161            0.109
#' #  Y3           0.779        0.180            0.110
#' #  Y4           0.728        0.039            0.097
#' #  Y5           0.283        0.109            0.431
#' #  Y6           0.412       -0.011            0.365
#' #  Y7           0.010        0.597           -0.150
#' #  Y8           0.516        0.009            0.195
#' #  Y9           0.578        0.092            0.057
#' #  Y10          0.010        0.820           -0.108
#' #  Y11         -0.043        0.696            0.119
#' #  Y12          0.024        0.652            0.222
#' #  Y13          0.376        0.123            0.350
#' #  Y14          0.073        0.194            0.639
#' #  Y15          0.052        0.183            0.693
#' #  Y16         -0.026        0.161            0.445
#' #  Y17          0.042        0.025            0.696
#' #  Y18         -0.019        0.175            0.607
#'
#'
#' @references Asparouhov, T., & Muthen, B. (2009). Exploratory structural equation modeling. Structural equation modeling: a multidisciplinary journal, 16(3), 397–438.
#' @references Mai, Y., Srivastava, D.K., & Krull, K.R. (2021). Estimating Composite reliability of Multidimensional Measurement with Overlapping Items. Present at the 2021 Eastern North American Region (ENAR) Spring Virtual Meeting.
#' @references Morin, A. J. S., Arens, A. K., & Marsh, H. W. (2016). A Bifactor Exploratory Structural Equation Modeling Framework for the Identification of Distinct Sources of Construct-Relevant Psychometric Multidimensionality. Structural equation modeling, 23(1), 116–139. doi: 10.1080/10705511.2014.961800
#' @references Varni, J. W., Burwinkle, T. M., Katz, E. R., Meeske, K., & Dickinson, P. (2002). The PedsQL in pediatric cancer: Reliability and validity of the Pediatric Quality of Life Inventory Generic Core Scales, Multidimensional Fatigue Scale, and Cancer Module. Cancer, 94(7), 2090.
#' @references St. Jude Children's Research Hospital. SJCRH. (2007-2021). St. Jude LIFE Study. \url{https://www.stjude.org/research/initiatives/cancer-survivorship-research/st-jude-life-study.html}
"PedsQLMFS"



