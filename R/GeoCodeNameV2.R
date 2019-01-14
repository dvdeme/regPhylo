#' @title Retrieve geographic coordinates from place names using nominatim openstreetmap API
#'
#' @description This function uses the nominatim openstreetmap API
#' (https://nominatim.openstreetmap.org) to retrieve
#' geographic coordinates for the sequences without geographic coordinates but
#' carrying information about their sampling location.
#'
#' @details Since the RGoogleMap API required a paying key (changed policy the July 16 2018)
#' we used the nominatim openstreetmap API.
#'
#' @param input a table coming from the function \code{\link{GeoCoord.WGS84}} or from
#'  \code{\link{Query.GeOMe.XY}}.
#' @param output name of the output table exported into the working directory.
#'
#' @param CorrTab an optional two column table providing a correction for the location names.
#' The first column provides the original location name, the second column provides a
#' corrected name fitting the purpose of the study.
#'
#' @param AutoCorrNZ if 'TRUE' then the correction of the place name will be applied
#' based on New Zealand (i.e. relevant for the phylogenetic tree of New Zealand marine
#' ray-finned fishes, Eme et al. Submitted). However other corrections will be more
#' suitable in other contexts.
#'
#' @return Two identical tables are exported, one in the R environment and one in the working directory.
#' These tables includes additional fields: Location_used= name
#' of the location used to infer the geographic coordinates,
#' Geo_accuracy = precision of the coordinates, with
#' 'Inferred' = imprecised geographic coordinates inferred from the Location name,
#' 'From_DB' = accurate geographic coordinates (provided by the databases), and
#' 'NoLocationFound' = the geographic coordinates cannot be inferred from the place
#' name.

#' @examples ## Import the example data coming from GeoCoord.WGS84 function.
#' \dontrun{
#' data(Seq.Diastocapen)
#' Seq.DF1 = Seq.Diastocapen$Seq.DF1
#'
#' # Run the function to retrieve geographic coordinates for sequences without coordinates.
#' Seq.DF3 = GeoCodeNameV2(input = Seq.DF1, output = "Seq.DF3.txt",
#' AutoCorrNZ = TRUE)
#'
#' # To see the modifications
#' Seq.DF3[, c(20,21, 23, 24, 25)]
#'
#' # Run the function including a correction table
#' # Build a dummy correction table, including the following correction.
#' correctionTab = matrix(c("Tasman Sea: , , NA, South Norfolk Ridge, off the Three Kings Islands",
#' "New Zealand: Three Kings Islands"), ncol=2)
#' colnames(correctionTab) = c("OriginalLocationName", "CorrectedLocationName")
#' head(correctionTab)
#'
#' # Run the function using the correction table built above, but without the automatic
#' # corrections set-up for the study of New Zealand marine fishes (i.e. AutoCorrNZ = FALSE)
#' Seq.DF3_test = GeoCodeNameV2(input = Seq.DF1, output = "Seq.DF3.txt",
#' CorrTab = correctionTab, AutoCorrNZ = FALSE)
#'
#' # To see the modifications
#' Seq.DF3_test[, c(20,21, 23, 24, 25)]
#'
#' }
#'
#' @export GeoCodeNameV2


GeoCodeNameV2 = function(input = NULL, output = NULL, CorrTab = NULL, AutoCorrNZ = NULL) {

  # the geocodeOpenStreetMap function allowing to query nominatim openstreetmap API.
  geocodeOpenStreetMap = function(input) {
    ### remove space , ";", ":", ",", ".", for the URLs
    input = gsub(" ", "%20", input) # space
    input = gsub(":", "%3A", input) # ":"
    input = gsub(";", "%3B", input) # ";"
    input = gsub(",", "%3C", input) # ","
    url  =paste("http://nominatim.openstreetmap.org/search?q=", input,
              "&limit=20&format=json", sep="") ### build complete url for the request
    # with maximun 20 answers.
    outres = RJSONIO::fromJSON(url)
    if(is.vector(outres)){
      res = c(outres[[1]]$lat, outres[[1]]$lon)
    } else {
      res = c(NA, NA)
    }
    return(res)
  }

  # First clean the comma and other empty cells in the field 'Location'.
    LocatCl = gsub(" :  ,  ,  ", NA, input[, 20], fixed = TRUE)
    LocatCl = gsub(" :  ,  [,]?[ ]?[ ]?", NA, LocatCl, fixed = TRUE)
    LocatCl = gsub("(,  ,[ ]?[ ]?[,]?$)", "", LocatCl, perl = TRUE)
    LocatCl = gsub(" :  ,  ,  ,[ ]+", NA, LocatCl, perl = TRUE)
    LocatCl = gsub("(,  ,)", ", ", LocatCl, perl = TRUE)

    # Keep track of the Original location name
    OrigLocatName = input[,20]

    # If a Correction table is proposed
    if(is.null(CorrTab) == FALSE){
      ### remove all the NA in the first column
      if(length(which(is.na(CorrTab[,1]) ==TRUE)) >0){
      CorrTab2 = CorrTab[-which(is.na(CorrTab[,1]) ==TRUE),]
      } else {
        CorrTab2 = CorrTab
      }
      i=1
      for(i in 1:dim(CorrTab2)[1]){
        LocatCl[which(LocatCl == CorrTab2[i,1])] = CorrTab2[i,2]
      }
    }

    input = cbind(input[, c(1:19)], Location = LocatCl, input[, c(21:26)], Oriorder = seq(1,
                                                                                          dim(input)[1], by = 1))

    # Select only the cells with information in the field Location but without the
    # precise geographic coordinates already recorded.
    pos1 = input[which(is.na(input[, 20]) == "FALSE"), ]
    pos.Seq1 = pos1[which(is.na(pos1[, 22]) == "TRUE"), 27]

    Locat = rep(NA, dim(input)[1])
    Locat[pos.Seq1] <- as.character(input[pos.Seq1, 20])

    # In the absence of geographic coordinates and Location, select the geographic
    # information brought by the field 'isolation_source'.
    pos.2 = input[which(is.na(input[, 21]) == "FALSE"), ]  ## Select only the rows without NA in the field 'isolation source'.
    pos.3 = pos.2[which(is.na(pos.2[, 22]) == "TRUE"), ]  ## Remove the rows with geographic coordinates.
    pos.4 = pos.3[which(is.na(pos.3[, 20]) == "TRUE"), 27]  ## Remove the rows with a location name in the field 'location'.

    Locat[pos.4] = as.character(input[pos.4, 21])  # Add the non NA isolation source to the Locat object.
    Locat = Locat[which(is.na(Locat) == "FALSE")]  # Remove the NA.
    pos.Seq = sort(c(pos.Seq1, pos.4))  # Keep the right order of the rows.

    # Clean the unnecessary NAs.
    Locat = gsub(", NA", "", Locat, fixed = TRUE)


    # Perform correction of the place names for the study focussing on the New Zealand marine ray-finned fishes.
    if(AutoCorrNZ == TRUE) {
    # Fastidious cleaning of the description of the locations tailored for our study.
    Locat = gsub("misc_structure [ ]?<1..>224", "", Locat, perl = TRUE)
    Locat = gsub("misc_feature", "", Locat, fixed = TRUE)
    Locat = gsub("([ ]+rRNA$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+tRNA$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+tRNA$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+CDS$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+mRNA$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+gene$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+D-loop$)", "", Locat, perl = TRUE)
    Locat = gsub("([ ]+misc_RNA$)", "", Locat, perl = TRUE)
    Locat = gsub("repeat_region", "", Locat, fixed = TRUE)
    Locat = gsub(":East Coast", "", Locat, fixed = TRUE)
    Locat = gsub(", 2184/120", "", Locat, fixed = TRUE)
    Locat = gsub("Slope water,", "", Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, Australian Antarctic Territory, Heard and Macdonald Islands.*",
        "Heard and Macdonald Islands", Locat, perl = TRUE)
    Locat = gsub("French Polynesia: French Polynesia, Society Islands.*", "French Polynesia: French Polynesia, Society Islands",
        Locat, perl = TRUE)
    Locat = gsub(".*Bay of Biscay", "Bay of Biscay", Locat, perl = TRUE)
    Locat = gsub("South Africa: South Africa, Agulhas Bank", "Agulhas Bank", Locat,
        fixed = TRUE)
    Locat = gsub("New Zealand: New Zealand, Chatham Rise.*", "New Zealand: Chatham Rise",
        Locat, perl = TRUE)
    Locat = gsub("New Zealand: New Zealand.*Norfolk Ridge.*", "New Zealand: Norfolk Ridge",
        Locat, perl = TRUE)
    Locat = gsub("Antarctica: Antarctica, South Sandwich Island/Scotia ridge complex.*",
        "South Sandwich Island", Locat, perl = TRUE)
    Locat = gsub("Antarctica: South Sandwich Island/Scotia ridge complex.*", "South Sandwich Island",
        Locat, perl = TRUE)
    Locat = gsub("New Zealand: New Zealand, Southern Plateau.*", "New Zealand", Locat,
        perl = TRUE)
    Locat = gsub("South Africa: South Africa.*west coast", "South Africa: South Africa, west coast",
        Locat, perl = TRUE)
    Locat = gsub("United States: United States, California, Eureka.*", "United States: California, Eureka",
        Locat, perl = TRUE)
    Locat = gsub("Antarctica: Antarctica, Ross Sea.*", "Ross Sea", Locat, perl = TRUE)
    Locat = gsub("France.*Western Mediterranean", "France: Collioure", Locat, perl = TRUE)
    Locat = gsub("Italy: Italy.*Adriatic Sea", "Italy: Adriatic Sea", Locat, perl = TRUE)
    Locat = gsub("Greece: Greece.*Aegean Sea.*", "Greece: Aegean Sea", Locat, perl = TRUE)
    Locat = gsub("United Kingdom: United Kingdom, Liverpool Bay.*", "United Kingdom: Liverpool Bay",
        Locat, perl = TRUE)
    Locat = gsub("New Zealand: New Zealand, TAN.*", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("Australia: Australia, Western Australia, Rottnest Island.*", "Australia: Western Australia, Rottnest Island",
        Locat, perl = TRUE)
    Locat = gsub("Australia: Australia.*Tasman Sea.*", "Australia: Tasman Sea", Locat,
        perl = TRUE)
    Locat = gsub("Atlantic Ocean: northern middle", "Atlantic Ocean: Acores", Locat,
        fixed = TRUE)
    Locat = gsub("Australia: Australia, Tasmania.*East of Maria Island", "Australia: Australia, Tasmania, East of Maria Island",
        Locat, perl = TRUE)
    Locat = gsub("United States: United States, Hawaii, Hawaii,Northwestern Hawaiian Islands, Maro Reef.*",
        "United States: Hawaii, Northwestern Hawaiian Islands, Maro Reef", Locat,
        perl = TRUE)
    Locat = gsub("Atlantic Ocean: Atlantic Ocean, Mid-Atlantic Bight", "Mid-Atlantic Bight",
        Locat, perl = TRUE)
    Locat = gsub("South Indian Sea", "", Locat, fixed = TRUE)
    Locat = gsub("Israel: Israel, Mediterranean", "Israel: Israel", Locat, fixed = TRUE)
    Locat = gsub("Chile: Chile, O'Higgins, Bajo O\\'Higgins", "Chile: Chile, Bajo O'Higgins",
        Locat, perl = TRUE)
    Locat = gsub("United States Virgin Islands: United States Virgin Islands, Saint Thomas.*",
        "United States: Virgin Islands, Saint Thomas", Locat, perl = TRUE)
    Locat = gsub("Atlantic Ocean: eastern", "Capo Verde", Locat, fixed = TRUE)
    Locat = gsub("USA: Hawaii, Northwestern Hawaiian Islands, Maro Reef.*", "USA: Hawaii, Northwestern Hawaiian Islands, Maro Reef",
        Locat, perl = TRUE)
    Locat = gsub(".*[Mm]id[ -]Atlantic [Bb]ight", "Mid Atlantic bight", Locat, perl = TRUE)
    Locat = gsub("USA.*St. Georges Bank.*", "USA: St. Georges Bank", Locat, perl = TRUE)
    Locat = gsub("New Zealand: Three Kings Island.*", "New Zealand: Three Kings Island",
        Locat, perl = TRUE)
    Locat = gsub("USA: Gulf of Mexico", "Arrecife, Alacranes island", Locat, fixed = TRUE)
    Locat = gsub("Turkey.*Iskenderun Bay", "Turkey: Iskenderun Bay", Locat, perl = TRUE)
    Locat = gsub("Southern Ocean.*South Tasman Rise.*", "South Tasman Rise", Locat,
        perl = TRUE)
    Locat = gsub("Pacific Ocean:Louisville Ridge.*", "Pacific Ocean:Louisville Ridge",
        Locat, perl = TRUE)
    Locat = gsub("Ireland: Rockall Trough", "Ireland", Locat, fixed = TRUE)
    Locat = gsub("Atlantic Ocean.*Atlantic Ocean", "Atlantic Ocean", Locat, perl = TRUE)
    Locat = gsub("Australia: Australia, Tasmania, Derwent Estuary", "Australia: Tasmania",
        Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, Western Australia, Point Samson.*", "Australia: Western Australia, Point Samson",
        Locat, perl = TRUE)
    Locat = gsub("Indian Ocean: Indian Ocean", "Indian Ocean", Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, Tasmania, Maatsuyker Hill", "Australia: Tasmania, Maatsuyker island",
        Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, Western Australia.*Shark Bay", "Australia: Western Australia, Shark Bay",
        Locat, perl = TRUE)
    Locat = gsub("Iran.*Bushehr, Nayband National Park Coast", "Iran: Bushehr, Nayband National Park Coast",
        Locat, perl = TRUE)
    Locat = gsub("South Africa: South Africa, KwaZulu-Natal.*", "South Africa: KwaZulu-Natal",
        Locat, perl = TRUE)
    Locat = gsub("Japan: Japan, Kanto, Yokosuka.*", "Japan: Kanto, Yokosuka", Locat,
        perl = TRUE)
    Locat = gsub("Australia: Australia, Tasmania, Macquarie Island.*", "Australia: Tasmania, Macquarie Island",
        Locat, perl = TRUE)
    Locat = gsub("Atlantic Ocean.*Flemish Cap", "Flemish Cap", Locat, perl = TRUE)
    Locat = gsub("South Africa: South Africa, Walters Shoal", "South Africa", Locat,
        fixed = TRUE)
    Locat = gsub("South Africa: South Africa, Umlaas Deep", "South Africa", Locat,
        fixed = TRUE)
    Locat = gsub("South Pacific Ocean: South Pacific Ocean, Southwest Pacific Ocean",
        "Pacific Ocean: Louisville Ridge", Locat, fixed = TRUE)
    Locat = gsub("Southern Ocean: Heard Island [S]?[E]?", "Heard Island", Locat,
        perl = TRUE)
    Locat = gsub("Japan:around Yaeyama Islands", "Yaeyama islands", Locat, fixed = TRUE)
    Locat = gsub("India: Arabian Sea", "Arabian Sea", Locat, fixed = TRUE)
    Locat = gsub("Denmark: North East Atlantic", "Faroe Islands", Locat, fixed = TRUE)
    Locat = gsub("Denmark: Denmark,  North East Atlantic,  ", "Faroe islands", Locat,
        fixed = TRUE)
    Locat = gsub("Atlantic Ocean: Mid-Atlantic [Rr]idge", "Azores island", Locat,
        perl = TRUE)
    Locat = gsub("Russia: Azov Sea, Tuzla Cape", "Sea of Azov ", Locat, fixed = TRUE)
    Locat = gsub("Russia: Sea of Japan, South Primorye, Vostok Bay", "Vostok Bay",
        Locat, fixed = TRUE)
    Locat = gsub("Italy: Mediterranean Sea, Sardinia, San Giovanni Lagoon", "Sardinia",
        Locat, fixed = TRUE)
    Locat = gsub("Ecuador: Bahia Divine,[ ]?Isla Santa Cruz, Galapagos Isld", "Ecuador: Galapagos islands, Isla Santa Cruz",
        Locat, perl = TRUE)
    Locat = gsub("Ecuador: Galapagos, Bahia Divine,[ ]?Isla Santa Cruz, Galapagos Isld",
        "Ecuador: Galapagos islands, Isla Santa Cruz", Locat, perl = TRUE)
    Locat = gsub("Iran:Caspian Sea", "Caspian Sea", Locat, fixed = TRUE)
    Locat = gsub("France: Clipperton Island, East Pacific", "Clipperton island",
        Locat, fixed = TRUE)
    Locat = gsub("Jordania:Red Sea", "Aqaba", Locat, fixed = TRUE)
    Locat = gsub("Jordan:Red Sea", "Aqaba", Locat, fixed = TRUE)
    Locat = gsub("Atlantic Ocean: equatorial", "Atlantic Ocean", Locat, fixed = TRUE)
    Locat = gsub("Tonga: Eua, Ha'Aluma Beach", "Tonga", Locat, fixed = TRUE)
    Locat = gsub("USA: Puerto Rico North Coast", "USA: Puerto Rico", Locat, fixed = TRUE)
    Locat = gsub("Cape Verde: Dacia-Concepcion", "Cape Verde", Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean: Rapa Nui", "Easter island", Locat, fixed = TRUE)
    Locat = gsub("USA: Raita, Hawaiian Archipelago", "Hawaii", Locat, fixed = TRUE)
    Locat = gsub("USA: Gardner, Hawaiian Archipelago", "Hawaii", Locat, fixed = TRUE)
    Locat = gsub("Australia: Seafood Trade, Sydney Fish.+", "Australia: Sydney",
        Locat, perl = TRUE)
    Locat = gsub("France: Reunion, Indian Ocean", "Reunion island", Locat, fixed = TRUE)
    Locat = gsub("Australia: Christmas Island[s]?[ ]?", "Christmas island", Locat,
        perl = TRUE)
    Locat = gsub("USA: French Frigate Shoals, HI", "French Frigate Shoals", Locat,
        fixed = TRUE)
    Locat = gsub("Antarctica: Antarctic Peninsula", "Antarctica: Graham land", Locat,
        fixed = TRUE)
    Locat = gsub("Antarctica: near Italy's Terra Nova Bay station", "Antarctica: Scott coast",
        Locat, fixed = TRUE)
    Locat = gsub("New Zealand: Karai Kari", "New Zealand: Karikari", Locat, fixed = TRUE)
    Locat = gsub("Australia: Sydney Harbor[ ]?", "Australia: Sydney", Locat, perl = TRUE)
    Locat = gsub("^.*Tambourine Bay.*$", "Australia: Tambourine Bay", Locat, perl = TRUE)
    Locat = gsub("Japan:Atlantic Ocean", "Japan", Locat, fixed = TRUE)
    Locat = gsub("Mauritania: Atlantic Ocean", "Mauritania", Locat, fixed = TRUE)
    Locat = gsub("Venezuela: Caribbean Sea", "Venezuela: Gulf of Venezuela", Locat,
        fixed = TRUE)
    Locat = gsub("Taiwan: southeastern coast", "Taiwan", Locat, fixed = TRUE)
    Locat = gsub("New Zealand: northeastern coast of North Island", "New Zealand: Poor Knights Island",
        Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean: Central Pacific", "Kiribati", Locat, fixed = TRUE)
    Locat = gsub("Taiwan: East China Sea, Off Suao", "Taiwan: Suao", Locat, fixed = TRUE)
    Locat = gsub("Taiwan: East China Sea, Off Keelun", "Taiwan: East China Sea, Off Keelun",
        Locat, fixed = TRUE)
    Locat = gsub("Israel: Israel,  , Levant Basin,", "Israel: Israel", Locat, fixed = TRUE)
    Locat = gsub("United States: United States, New York, New York City.*$", "United States: New York, New York City",
        Locat, perl = TRUE)
    Locat = gsub("Antarctica: Antarctica,  ,  , Burdwood Banks", "Burdwood Banks",
        Locat, fixed = TRUE)
    Locat = gsub("United States: United States, New Jersey, Barnegat Light, Cassidy's Fish Market",
        "United States: New Jersey, Barnegat Light", Locat, fixed = TRUE)
    Locat = gsub("India: India, Tamil Nadu, Rameshwaram, Kilakarai coast", "India: Tamil Nadu, Rameshwaram",
        Locat, fixed = TRUE)
    Locat = gsub("Cyprus: Cyprus,  , Levant Basin,  ", "Cyprus", Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean: Pacific Ocean,  , Central Eastern Pacific,  ", "Clipperton island",
        Locat, fixed = TRUE)
    Locat = gsub(" :  ,  , E-02-04_L.70", NA, Locat, fixed = TRUE)
    Locat = gsub("Mexico: Mexico,  ,  , Hurricane Bank; from stomach of wahoo", "Mexico: Cabo San Lucas",
        Locat, fixed = TRUE)
    Locat = gsub("New Zealand TRIP2626/118", "New Zealand", Locat, fixed = TRUE)
    Locat = gsub("(T[Aa][Nn]0.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Chatham Rise, Chatham Rise.*$)", "New Zealand", Locat,
        perl = TRUE)
    Locat = gsub("(OBS[ ]?2.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("(KAH0.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Challenger Plateau.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("New Zealand, OBS 227/1002, Challenger Plateau", "New Zealand",
        Locat, fixed = TRUE)
    Locat = gsub("New Zealand, Ocean survey 20/20 TAN 0707/056: RV Tangaroa: Challenger Plateau",
        "New Zealand", Locat, fixed = TRUE)
    Locat = gsub("(New Zealand, Wanganella Bank.*$)", "Norfolk Island", Locat, perl = TRUE)
    Locat = gsub("New Zealand, Marlborough Sounds, Marlborough Sounds", "New Zealand: Marlborough Sounds",
        Locat, fixed = TRUE)
    Locat = gsub("(SWA0.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Ross Sea, Ross Sea.*$)", "Ross Sea", Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Lord Howe Rise.*$)", "Lord Howe island", Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Chatham Rise.*$)", "Chatham island", Locat, perl = TRUE)
    Locat = gsub("New Zealand, South-east Chatham Slope area.*$", "Chatham island",
        Locat, perl = TRUE)
    Locat = gsub("(New Zealand, Just North of the Moko Hinau Islands.*$)", "New Zealand, Moko Hinau Islands",
        Locat, perl = TRUE)
    Locat = gsub("New Zealand, TAN, TAN,", "New Zealand", Locat, fixed = TRUE)
    Locat = gsub("(New Zealand, ECSI, ECSI.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("New Zealand, Trip29.*$", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("(New Zealand THH05.*$)", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("^.*Techobanine$", "South Africa", Locat, perl = TRUE)
    Locat = gsub("(United States: United States, Hawaii, Hawaii,Northwestern Hawaiian Islands.*&)",
        "Hawaii islands", Locat, perl = TRUE)
    Locat = gsub("Denmark: Denmark,  , North East Atlantic, ", "Faroe islands", Locat,
        fixed = TRUE)
    Locat = gsub("(Antarctica: Antarctica,  , South Sandwich Island/Scotia.*$)",
        "Sandwich Island", Locat, perl = TRUE)
    Locat = gsub("New Zealand: New Zealand,  , Campbell Plateau, ", "New Zealand",
        Locat, perl = TRUE)
    Locat = gsub("Russia: Russia,  , Azov Sea,  ", "Azov Sea", Locat, fixed = TRUE)
    Locat = gsub("(Ecuador: Ecuador, Galapagos.*$)", "Galapagos islands", Locat,
        perl = TRUE)
    Locat = gsub("Australia: Australia, New South Wales, Seafood Trade, Sydney Fish Markets, NSW,",
        "Australia: Australia, Sydney", Locat, fixed = TRUE)
    Locat = gsub("(Italy: Italy,  [,]?[ ]?Cabras Lagoon.*$)", "Italy: Sardinia, Oristano",
        Locat, perl = TRUE)
    Locat = gsub("South Africa: South Africa, Mid-Indian Ocean", "Indian Ocean",
        Locat, fixed = TRUE)
    Locat = gsub("South Africa: South Africa,  ,  , Walters Shoal", "South Africa",
        Locat, fixed = TRUE)
    Locat = gsub("South Africa: South Africa, NA,  Walters Shoal", "South Africa",
        Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, Western Australia, somewhere in the Western Deepwater Trawl,",
        "Australia: Australia, Western Australia", Locat, fixed = TRUE)
    Locat = gsub(" :  ,  ,  , Aquarium trade specimens, capture locality unknown",
        "", Locat, fixed = TRUE)
    Locat = gsub("United States: United States, Hawaii, Hawaii,Oahu, Kailua Beach, stranded on shore,",
        "Hawaii islands", Locat, fixed = TRUE)
    Locat = gsub("Portugal: Portugal,  [,]?[ ]?Praia do Baleal,", "Portugal, Baleal",
        Locat, perl = TRUE)
    Locat = gsub("Departement code 75", "Paris", Locat, fixed = TRUE)
    Locat = gsub("Departement code 13", "Marseille", Locat, fixed = TRUE)
    Locat = gsub("Departement code 86", "Poitiers", Locat, fixed = TRUE)
    Locat = gsub("Atlantic Ocean: FAO 27, Galician Bank", "Spain: Galicia", Locat,
        fixed = TRUE)
    Locat = gsub("Pacific Ocean:western North Pacific Ocean", "Kuril islands", Locat,
        fixed = TRUE)
    Locat = gsub("Atlantic Ocean: mid-Atlantic, ridge", "Azores", Locat, fixed = TRUE)
    Locat = gsub("Gulf of Biscay (Northern Atlantic ocean)", "Bay of Biscay", Locat,
        fixed = TRUE)
    Locat = gsub("^Atlantic Ocean: south[ ]?$", "Bouvet island", Locat, perl = TRUE)
    Locat = gsub("USA: O'ahu, Hawaiian islands", "O'ahu island", Locat, fixed = TRUE)
    Locat = gsub("Australia: Cape Naturaliste.*$", "Australia: Cape Naturaliste",
        Locat, perl = TRUE)
    Locat = gsub("Antarctica: near Italy's Terra Nova Bay station", "Antarctica",
        Locat, fixed = TRUE)
    Locat = gsub("Atlantic Ocean: Eastern Central", "Capo Verde", Locat, fixed = TRUE)
    Locat = gsub("Japan: Nishikihoku", "Japan", Locat, fixed = TRUE)
    Locat = gsub("Japan:Wakayama, Nachi", "Japan: Wakayama", Locat, fixed = TRUE)
    Locat = gsub("Venezuela: Cayo de Agua, Los Roques islands", "Los Roques islands",
        Locat, fixed = TRUE)
    Locat = gsub("Taiwan: East China Sea, Off Keelun", "Taiwan: Keelun", Locat, fixed = TRUE)
    Locat = gsub("Israel: Israel,  Levant Basin.*$", "Israel", Locat, perl = TRUE)
    Locat = gsub("Antarctica: Antarctica,   , Burdwood Banks", "Antarctica", Locat,
        fixed = TRUE)
    Locat = gsub("Antarctica: Antarctica,   , Kerguelen Isl.", "Kerguelen island",
        Locat, fixed = TRUE)
    Locat = gsub("Cyprus: Cyprus,  Levant Basin, ", "Cyprus", Locat, fixed = TRUE)
    Locat = gsub("New Zealand, OBS.*$", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("New Zealand, WCSI, WCSI.*$", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("New Zealand, ton0.*$", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("New Zealand, Trip2.*$", "New Zealand", Locat, perl = TRUE)
    Locat = gsub("United States: United States, Hawaii, Hawaii,Northwestern Hawaiian Islands, Gardner Pinnacles, 30 m depth,  ",
        "Hawaii", Locat, fixed = TRUE)
    Locat = gsub("Turkey: Turkey,  NA", "Turkey", Locat, fixed = TRUE)
    Locat = gsub("New Zealand: New Zealand,  Southern Plateau, New Zealand", "New Zealand: Stewart island",
        Locat, fixed = TRUE)
    Locat = gsub("Antarctica: Antarctica,  South Sandwich Island/Scotia ridge complex, .*$",
        "Sandwich Island", Locat, perl = TRUE)
    Locat = gsub("New Zealand: New Zealand,  Campbell Plateau,  ", "New Zealand",
        Locat, fixed = TRUE)
    Locat = gsub("Australia: Australia, New South Wales, Seafood Trade, Sydney Fish Markets, NSW",
        "Sydney", Locat, fixed = TRUE)
    Locat = gsub("South Pacific Ocean: South Pacific Ocean,  Southwest Pacific Ocean, NA",
        "Chatham island", Locat, fixed = TRUE)
    Locat = gsub("North Pacific Ocean: North Pacific Ocean, NA, Northwest Pacific Ocean, NA",
        "Kuril islands", Locat, fixed = TRUE)
    Locat = gsub("French Polynesia: French Polynesia,  Society Islands, Haapiti Pass",
        "Society Islands", Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean: Pacific Ocean,  Central,  ", "Kiribati", Locat,
        fixed = TRUE)
    Locat = gsub("Pacific Ocean: Pacific Ocean,  Western,  ", "Marshall islands",
        Locat, fixed = TRUE)
    Locat = gsub("South Korea: South Korea, Jeju, Jeju Island", "Jeju island", Locat,
        fixed = TRUE)
    Locat = gsub("^Pacific Ocean$", "Kiribati", Locat, perl = TRUE)
    Locat = gsub("Pacific Ocean: Pacific Ocean,  Off the coast of New Caledonia",
        "New Caledonia", Locat, fixed = TRUE)
    Locat = gsub("Taiwan: Taiwan,  NA", "Taiwan", Locat, fixed = TRUE)
    Locat = gsub(" :  ,  ,  ", NA, Locat, fixed = TRUE)
    Locat = gsub(" :  ", NA, Locat, fixed = TRUE)
    Locat = gsub(" :  ,  [,]?[ ]?[ ]?", NA, Locat, perl = TRUE)
    Locat = gsub("(,  ,[ ]?[ ]?[,]?$)", "", Locat, perl = TRUE)
    Locat = gsub(" :  ,  ,  ,[ ]+", NA, Locat, perl = TRUE)
    Locat = gsub("(,  ,)", ", ", Locat, perl = TRUE)
    Locat = gsub("Japan:Kagoshima, Yaku", "Japan: Kagoshima", Locat, fixed = TRUE)
    Locat = gsub("USA: Bering Sea, Alaska", "Bering Sea", Locat, fixed = TRUE)
    Locat = gsub("France: Marquesas Islands", "Marquesas Islands", Locat, fixed = TRUE)
    Locat = gsub("USA: Newport River and Adams Creek, North Carolina", "USA: North Carolina, Newport River",
        Locat, fixed = TRUE)
    Locat = gsub("USA: Chocolate Bay and West Bay in Western Galveston Bay, Texas",
        "USA: Texas, Galveston Bay", Locat, fixed = TRUE)
    Locat = gsub("USA: Kailua Bay, Oahu, Hawaii ", "USA: Hawaii, Oahu", Locat, fixed = TRUE)
    Locat = gsub("France: Collioure Sea, Gulf of Lyon, Sete", "France: Sete", Locat,
        fixed = TRUE)
    Locat = gsub("Japan:Kanagawa, off Jogashima Island", "Japan: Jogashima Island",
        Locat, fixed = TRUE)
    Locat = gsub("Ecuador: Galapagos, Bahia Divine, Isla Santa Cruz, .*$", "Ecuador: Galapagos islands, Isla Santa Cruz",
        Locat, perl = TRUE)
    Locat = gsub("[ ]+marine waters", "", Locat, perl = TRUE)
    Locat = gsub("stomach of Nototodarus gouldi", NA, Locat, fixed = TRUE)
    Locat = gsub("GAB SARDI Survey Station MI2", NA, Locat, fixed = TRUE)
    Locat = gsub("Eumetopias jubatus (Stellar sealion) soft scat", NA, Locat, fixed = TRUE)
    Locat = gsub("stomach contents", NA, Locat, fixed = TRUE)
    Locat = gsub("^fish market$", NA, Locat, perl = TRUE)
    Locat = gsub("blind sample [HDS]C for ring trial", NA, Locat, perl = TRUE)
    Locat = gsub("^leaf$", NA, Locat, perl = TRUE)
    Locat = gsub("^seafood trade$", NA, Locat, perl = TRUE)
    Locat = gsub("individual from the ", "", Locat, fixed = TRUE)
    Locat = gsub("CSIRO, Hobart unregistered specimen trawled off eastern Tasmania by G. Yearsley",
        "Tasmania: Bicheno", Locat, fixed = TRUE)
    Locat = gsub("[ ]+by G. Yearsley", "", Locat, fixed = TRUE)
    Locat = gsub("[Aa]rchipelago", "islands", Locat, perl = TRUE)
    Locat = gsub("Galicia Bank [(]?Northern Atlantic.*$", "Spain: Galicia", Locat,
        perl = TRUE)
    Locat = gsub("Banc Lansdowne, New Caledonia.*$", "New Caledonia", Locat, perl = TRUE)
    Locat = gsub(".*off Recif Jouan, Lifou Island.*$", "Lifou island", Locat, perl = TRUE)
    Locat = gsub(".*Southeast Victoria by G.K. Yearsley.*$", "Australia: Victoria",
        Locat, perl = TRUE)
    Locat = gsub("Western Mediterranean", "Balearic islands", Locat, perl = TRUE)
    Locat = gsub("Swedish exclusive economic zone", "Sweden", Locat, fixed = TRUE)
    Locat = gsub("Indo-Pacific [Oo]cean", "Sulawesi", Locat, perl = TRUE)
    Locat = gsub("Indo-Pacific [Ww]aters", "Sulawesi", Locat, perl = TRUE)
    Locat = gsub("Parangipettai coast", "India: Parangipettai", Locat, fixed = TRUE)
    Locat = gsub("aquacultured", NA, Locat, fixed = TRUE)
    Locat = gsub("GAB SARDI Survey Station B4", NA, Locat, fixed = TRUE)
    Locat = gsub("bloodfeeding parasites (Gnathia sp.)", NA, Locat, fixed = TRUE)
    Locat = gsub("^coral reef$", NA, Locat, perl = TRUE)
    Locat = gsub("Chinese fish market", "China", Locat, fixed = TRUE)
    Locat = gsub("^gut content$", NA, Locat, perl = TRUE)
    Locat = gsub("^ocean water$", NA, Locat, perl = TRUE)
    Locat = gsub("commercial food product", NA, Locat, perl = TRUE)
    Locat = gsub("^Atlantic$", "Atlantic Ocean", Locat, perl = TRUE)
    Locat = gsub("^Mediterranean$", "Mediterranean Sea", Locat, perl = TRUE)
    Locat = gsub("Mauritania coast", "Nouakchott", Locat, perl = TRUE)
    Locat = gsub("Canaries_Saharian Bank", "Canaries islands", Locat, perl = TRUE)
    Locat = gsub("^imported specimen purchased in a fish market$", NA, Locat, perl = TRUE)
    Locat = gsub("^Pacific$", "Pacific Ocean", Locat, perl = TRUE)
    Locat = gsub("^western South Atlantic Ocean$", "Tristan Da Cunha island", Locat,
        perl = TRUE)
    Locat = gsub("Gulf St Vincent", "Australia: Adelaide", Locat, perl = TRUE)
    Locat = gsub("Arabian Sea, 100m depth", "Arabian Sea", Locat, perl = TRUE)
    Locat = gsub("^coral reef-associated$", NA, Locat, perl = TRUE)
    Locat = gsub("Falkland Islands, Islas Malvinas", "Falkland Islands", Locat, fixed = TRUE)
    Locat = gsub("Antarctic region: Elephant Island", "Elephant island", Locat, fixed = TRUE)
    Locat = gsub("Elephant Island", "Elephant island", Locat, fixed = TRUE)
    Locat = gsub("^Cantabric Sea$", "Bay of Biscay", Locat, perl = TRUE)
    Locat = gsub("formalin treated sample", NA, Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean: Pacific Ocean,  Central Eastern Pacific,  ", "Clipperton island",
        Locat, fixed = TRUE)
    Locat = gsub(" ,", "", Locat, fixed = TRUE)
    Locat = gsub(", $", "", Locat, perl = TRUE)
    Locat = gsub(": [no locality data]", NA, Locat, fixed = TRUE)
    Locat = gsub("misc_structure <1..>224", "", Locat, fixed = TRUE)
    Locat = gsub("Italy:Sardinia-Loc.Tortoli", "Italy: Sardinia, Tortoli", Locat,
        fixed = TRUE)
    Locat = gsub("Kiribati: Kiritimati .*$", "Kiribati: Kiritimati", Locat, perl = TRUE)
    Locat = gsub(", Region [0-9]", "", Locat, perl = TRUE)
    Locat = gsub(", E. Point", "", Locat, fixed = TRUE)
    Locat = gsub(", IV-A", "", Locat, fixed = TRUE)
    Locat = gsub(", Middle Gable, NZ", "", Locat, fixed = TRUE)
    Locat = gsub(", HI", "", Locat, fixed = TRUE)
    Locat = gsub(", IV-B", "", Locat, fixed = TRUE)
    Locat = gsub("Luzon/", "", Locat, fixed = TRUE)
    Locat = gsub("^: $", NA, Locat, perl = TRUE)
    Locat = gsub("Peru: 18,5\032\032-20\032\032S - 80\032\032-90\032\032W, Triangle South", "Peru", Locat,
        fixed = TRUE)
    Locat = gsub("United States: New York, New York City,.*$", "United States: New York, New York City",
        Locat, perl = TRUE)
    Locat = gsub("from stomach of wahoo", "", Locat, fixed = TRUE)
    Locat = gsub(" TAN .*$", "", Locat, perl = TRUE)
    Locat = gsub(" Trip 2626/149", NA, Locat, fixed = TRUE)
    Locat = gsub(", Halipro2, 95", "", Locat, fixed = TRUE)
    Locat = gsub(", 30 m depth", "", Locat, fixed = TRUE)
    Locat = gsub("FAO 27-IXa", "", Locat, fixed = TRUE)
    Locat = gsub(", AVRO CHIEF 39", "", Locat, fixed = TRUE)
    Locat = gsub(", Sub-area IXa", "", Locat, fixed = TRUE)
    Locat = gsub(" NA, ", "", Locat, fixed = TRUE)
    Locat = gsub(": NA$", "", Locat, fixed = TRUE)
    Locat = gsub("Argentina: Southern Patagonia, off Santa Cruz", "Argentina: Southern Patagonia, Santa Cruz",
        Locat, fixed = TRUE)
    Locat = gsub(", J.L. Lenteng Agung No.9", "", Locat, fixed = TRUE)
    Locat = gsub(", SE side, 55 fathoms depth", "", Locat, fixed = TRUE)
    Locat = gsub(", somewhere in the Western Deepwater Trawl", "", Locat, fixed = TRUE)
    Locat = gsub(": Aquarium trade specimens, capture locality unknown", "", Locat,
        fixed = TRUE)
    Locat = gsub(", 99-899 Iwaena St. No. 103", "", Locat, fixed = TRUE)
    Locat = gsub("Departement code 44", "Loire-Atlantique", Locat, fixed = TRUE)
    Locat = gsub("Departement code 94", "Val-de-Marne", Locat, fixed = TRUE)
    Locat = gsub("Departement code 33", "Bordeaux", Locat, fixed = TRUE)
    Locat = gsub("Departement code 50", "Saint-Lo", Locat, fixed = TRUE)
    Locat = gsub("Departement code 77", "Seine-et-Marne", Locat, fixed = TRUE)
    Locat = gsub("Departement code 91", "Essonne", Locat, fixed = TRUE)
    Locat = gsub(" (East Coast)", "", Locat, fixed = TRUE)
    Locat = gsub("Mexico: 8 mile radius of Cabo San Lucas,.*$", "Mexico: Cabo San Lucas",
        Locat, perl = TRUE)
    Locat = gsub("Mexico: Cabo San Lucas, Mexico", "Mexico: Cabo San Lucas", Locat,
        fixed = TRUE)
    Locat = gsub("stranded on shore", "", Locat, fixed = TRUE)
    Locat = gsub("E51 Nouyen Oanh St.", "", Locat, fixed = TRUE)
    Locat = gsub(" off Futto", "", Locat, fixed = TRUE)
    Locat = gsub(" [Oo]ff ", "", Locat, perl = TRUE)
    Locat = gsub(":[Oo]ff ", ": ", Locat, perl = TRUE)
    Locat = gsub("Aegean Sea Bali", "Aegean Sea", Locat, fixed = TRUE)
    Locat = gsub("Cat Bi Fish Markets", "", Locat, fixed = TRUE)
    Locat = gsub("Open Sea", "", Locat, fixed = TRUE)
    Locat = gsub("Supermarket", "", Locat, fixed = TRUE)
    Locat = gsub("2397/25", "", Locat, fixed = TRUE)
    Locat = gsub("2844/9", "", Locat, fixed = TRUE)
    Locat = gsub(" $", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Saitama, Toda, Arakawa River", "Japan:Saitama, Toda", Locat,
        fixed = TRUE)
    Locat = gsub("Levant Basin", "", Locat, fixed = TRUE)
    Locat = gsub("Spain: Azores Islands", "Spain: Azores", Locat, fixed = TRUE)
    Locat = gsub("Ahuriri Point", "", Locat, fixed = TRUE)
    Locat = gsub("Tasmania:Horseshoe Island", "Tasmania", Locat, fixed = TRUE)
    Locat = gsub("Horseshoe Island, Tasmania", "Tasmania", Locat, fixed = TRUE)
    Locat = gsub("Japan:Aomori,.*", "Japan:Aomori", Locat, perl = TRUE)
    Locat = gsub("Japan:Wakayama,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Kyoto,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Iwate,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Kochi,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Chiba,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Japan:Miyazaki,.*$", "", Locat, perl = TRUE)
    Locat = gsub("Panama: Fort Randolph", "Panama: Colon", Locat, fixed = TRUE)
    Locat = gsub("Southern, Mediterranean,", "", Locat, fixed = TRUE)
    Locat = gsub("Antarctica: Burdwood Banks", "Burdwood Banks", Locat, fixed = TRUE)
    Locat = gsub("Coastal area, Kilakarai coast", "", Locat, fixed = TRUE)
    Locat = gsub("West Coast[,]? West Coast Southern Island", "", Locat, perl = TRUE)
    Locat = gsub("West Coast Southern Island, Okarito lagoon, DOC", "", Locat, fixed = TRUE)
    Locat = gsub("Southern Plateau, ", "", Locat, fixed = TRUE)
    Locat = gsub("Tohoku-chiho", "Tohoku", Locat, fixed = TRUE)
    Locat = gsub(": NA$", "", Locat, perl = TRUE)
    Locat = gsub(" ,", "", Locat, fixed = TRUE)
    Locat = gsub("^$", NA, Locat, perl = TRUE)
    Locat = gsub("Greece: Vontas", "Greece", Locat, fixed = TRUE)
    Locat = gsub("Vetch's Pier", "", Locat, fixed = TRUE)
    Locat = gsub("Japan: Iriomote Island, Kuira River", "Japan: Iriomote Island",
        Locat, fixed = TRUE)
    Locat = gsub("Tasman Sea: South Norfolk Ridge,the Three Kings Islands", "New Zealand: Three Kings Islands", Locat, fixed=TRUE)
    Locat = gsub("Gulf of Iskenderun", "Dortyol", Locat, fixed = TRUE)
    Locat = gsub("Italy: Cabras Lagoon", "Italy: Sardinia, Oristano", Locat, fixed = TRUE)
    Locat = gsub("Antarctica: Scott coast", "Ross Sea", Locat, fixed = TRUE)
    Locat = gsub("JapanIburi, Kochi", "Japan: Iburi, Kochi", Locat, fixed = TRUE)
    Locat = gsub("Antarctica: South Sandwich Island/Scotia ridge complex.*", "South Sandwich Island",
        Locat, perl = TRUE)
    Locat = gsub("French Polynesia: Society Islands, Moorea, Haapiti Pass", "French Polynesia: Moorea",
        Locat, fixed = TRUE)
    Locat = gsub("Pacific Ocean:the coast of New Caledonia", "New Caledonia", Locat,
        fixed = TRUE)
    Locat = gsub("commercially purchased caviar eggs", NA, Locat, fixed = TRUE)
    Locat = gsub("Rib Reef, GBR", "Great Barrier Reef", Locat, fixed = TRUE)
    Locat = gsub("South Atlantic Ridge hydrothermal vents", "South Atlantic", Locat, fixed = TRUE)
    Locat = gsub("Antarctic region: near New Zealand", "New Zealand: Campbell Island", Locat, fixed = TRUE)
    Locat = gsub("Zanzibar, Indian Ocean", "Zanzibar", Locat, fixed = TRUE)
    Locat = gsub("Tropical North Pacific", "USA: Hawaii", Locat, fixed = TRUE)
    Locat = gsub("western North Atlantic Ocean", "North Atlantic Ocean", Locat, fixed = TRUE)

    # Perfomed generic correction when the AutoCorrNZ is disabled.
    } else {
      Locat = gsub("misc_structure[ ]?[ ]?<1..>224", "", Locat, perl = TRUE)
      Locat = gsub("misc_feature", "", Locat, fixed = TRUE)
      Locat = gsub("([ ]+rRNA$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+tRNA$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+tRNA$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+CDS$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+mRNA$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+gene$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+D-loop$)", "", Locat, perl = TRUE)
      Locat = gsub("([ ]+misc_RNA$)", "", Locat, perl = TRUE)
      Locat = gsub("repeat_region", "", Locat, fixed = TRUE)
      Locat = gsub(" :  ,  ,  ", "", Locat, fixed = TRUE)
      Locat = gsub(" :  ,  [,]?[ ]?[ ]?", "", Locat, perl = TRUE)
      Locat = gsub("(,  ,[ ]?[ ]?[,]?$)", "", Locat, perl = TRUE)
      Locat = gsub(" :  ,  ,  ,[ ]+", "", Locat, perl = TRUE)
      Locat = gsub("(,  ,)", ", ", Locat, perl = TRUE)
      Locat = gsub("^: $", "", Locat, perl = TRUE)
      Locat = gsub(" ,", "", Locat, fixed = TRUE)
      Locat = gsub(", $", "", Locat, perl = TRUE)
      Locat = gsub(": [no locality data]", "", Locat, fixed = TRUE)
      Locat = gsub(" $", "", Locat, perl = TRUE)
      Locat = gsub(": NA$", "", Locat, perl = TRUE)
    }

    LocatTab = cbind(Locat, orderLocat = seq(1, length(Locat), by = 1))

    LocatUniq = unique(Locat)

    Latitude = vector()
    Longitude = vector()
    i = 1
    for (i in 1:length(LocatUniq)) {
        if (is.na(LocatUniq[i])) {
            a = c(NA, NA)
        } else {
          a = geocodeOpenStreetMap(LocatUniq[i])
        }
        Latitude = c(Latitude, a[[1]])
        Longitude = c(Longitude, a[[2]])
    }  # End for i.

    LocatDF = cbind(LocatUniq, Latitude, Longitude)
    LocatGeo = merge(LocatTab, LocatDF, by.x = 1, by.y = 1, all.x = T)
    LocatGeoOrd = LocatGeo[order(as.numeric(as.character(LocatGeo[, 2]))), ]  # Retrieve the initial order

    # Create a vector reporting the names of the location that have been used to
    # retrieve the geographic coordinates.
    Location_used = rep(NA, dim(input)[1])
    Location_used[pos.Seq] = as.character(LocatGeoOrd[, 1])

    # The geographic coordinates.
    Latitude_Y = input[, 22]
    Longitude_X = input[, 23]
    Latitude_Y[pos.Seq] = as.numeric(as.character(LocatGeoOrd[, 3]))
    Longitude_X[pos.Seq] = as.numeric(as.character(LocatGeoOrd[, 4]))

    # Precision of the coordinates, 'Inferred'= not precise inferred from the
    # Location name, 'From_DB'= more accurate.
    Geo_accuracy = rep(NA, dim(input)[1])
    Geo_accuracy[pos.Seq[which(is.na(Locat) == FALSE)]] = "Inferred"
    Geo_accuracy[which(is.na(input[, 22]) == "FALSE")] = "From_DB"

    # Replacing 'Inferred' from the column 'Geo_accuracy' by 'NoLocationFound', when
    # the location couldn't be found and the geographic coordinates are NA.
    Geo_accuracy[grep("Inferred", Geo_accuracy)[which(is.na(Latitude_Y[grep("Inferred",
        Geo_accuracy)]) == TRUE)]] = "NoLocationFound"

    # The output Table
    outputDF = cbind(input[, c(1:19)], Location = OrigLocatName, Location_used, isolation_source = input[,
        21], Latitude_Y, Longitude_X, Geo_accuracy, input[, c(24:26)])
    utils::write.table(outputDF, file = output, sep = "\t", row.names = FALSE)
    return(outputDF)
}  # End of the function.
