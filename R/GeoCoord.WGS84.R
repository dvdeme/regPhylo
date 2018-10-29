#' @title Extract longitude and latitude coordinates in separate fields and conversion in WGS84 using +- North, and +- East

#' @description This function splits the field 'Lat_Lon' in two distinct fields for 'Latitude'
#' and 'Longitude' and transform them to WGS84 geographic coordinates using +- North, and +- East (i.e. -45.10
#' -53.23)

#' @param input table coming from the function Congr.NCBI.BOLD.perReposit()
#' @param output name of the output table exported in the working directory.
#'
#' @details This function uses some code from the '05.extract_coords_from_GB.r' available at
#' 'https://github.com/paolo-gratton/Gratton_et_al_JBiogeogr_2016/blob/master/05.extract_coords_from_GB.r'
#' and used in the paper Gratton et al. 2016 'A world of sequences: can we use
#' georeferenced nucleotide databases for a robust automated phylogeography?'
#' Journal of Biogeography. DOI: 10.1111/jbi.12786.
#'
#' @return Two identical tables are exported, one in the R environment and one in the working directory.
#' These tables are similar to the input table except that they have one additional column,
#' the previous "Lat_Lon" corresponds now to two distinct columns "Latitude" and "Longitude".

#' @examples # Import a tables exported by the function Congr.NCBI.BOLD.perReposit
#' data(Seq.Diastocapen)
#' Seq.DF = Seq.Diastocapen$Seq.DF
#'
#' # Run the function
#' Seq.DF1=GeoCoord.WGS84(input=Seq.DF, output="Seq.DF.txt")
#' dim(Seq.DF1)

#' @export GeoCoord.WGS84


GeoCoord.WGS84 = function(input = NULL, output = NULL) {

    # Remove the NA NA in the dataset introduced when retrieving the BOLD
    # information.
    coord = gsub("NA NA", NA, input[, 22], fixed = TRUE)

    # Create two empty vectors for the new geographic coordinates in decimal degrees
    # (WGS84).
    all_lats = rep(NA, length(coord))
    all_lons = rep(NA, length(coord))

    # Regular expression for DECIMAL DEGREES with NSEW (with or without comma...)
    decimal_pattern_1 = "[0-9]{1,3}\\.[0-9]{0,18}\\s*?[NS],*?\\s[0-9]{1,3}\\.[0-9]{0,18}\\s*?[EW]"
    length(grep(decimal_pattern_1, coord))/length(input[, 22])
    # Regular expression for DECIMAL DEGREES with NSEW (without decimal= no point,
    # example 45 N 4 E)
    decimal_pattern_2 = "[0-9]{1,3}\\s*?[NS],*?\\s[0-9]{1,3}\\s*?[EW]"
    NSEW.pos = sort(c(grep(decimal_pattern_1, coord), grep(decimal_pattern_2, coord)))
    if (length(NSEW.pos) > 0)
        {
            # Check if the pattern is present in the datafile.
            coordNSEW = coord[NSEW.pos]
            lat_sign = rep("+", length(coordNSEW))
            lon_sign = rep("+", length(coordNSEW))
            lat_sign[grep("S", coordNSEW)] = "-"
            lon_sign[grep("W", coordNSEW)] = "-"
            ll_ll = strsplit(coordNSEW, "\\sN|\\sS|\\sE|\\sW")
            lat = as.numeric(sapply(ll_ll, "[[", 1))
            lon = as.numeric(sapply(ll_ll, "[[", 2))
            lat[lat_sign == "-"] = -lat[lat_sign == "-"]
            lon[lon_sign == "-"] = -lon[lon_sign == "-"]
            all_lats[NSEW.pos] = lat
            all_lons[NSEW.pos] = lon
        }  # End if(length(NSEW.pos)>0).


    # Regular expression for DECIMAL DEGREES without NSEW
    decimal_pattern_3 = "-*?[0-9]{1,3}\\.?[0-9]{0,10}[ ,]\\s*?-*?[0-9]{1,3}\\.?[0-9]{0,10}"
    # Select only the geographic coordinates following the pattern 3 = Decimal
    # coordinate with +-.
    coordDD = coord[grep(decimal_pattern_3, coord)]
    if (length(coordDD) > 0)
        {
            # Check if the pattern is present in the datafile.
            ll_ll = strsplit(coordDD, ",|, | ")
            lat = as.numeric(sapply(ll_ll, "[[", 1))
            lon = as.numeric(sapply(ll_ll, "[[", 2))
            all_lats[grep(decimal_pattern_3, coord)] = lat
            all_lons[grep(decimal_pattern_3, coord)] = lon
        }  # End if(length(coordDD)>0).


    # Write lat and lon for all coordinates with degrees and DECIMAL MINUTES.
    # Regular expression for degrees and DECIMAL MINUTES
    degmin_pattern_1 = "[0-9]{1,3}\\s*?([Dd]eg|[Dd]egrees)\\s*?[0-9]{1,2}.[0-9]{1,6}'\\s*?[NS][,;]\\s*?[0-9]{1,3}\\s*?([Dd]eg|[Dd]egrees)\\s*?[0-9]{1,2}.[0-9]{1,6}'\\s*?[EW]"
    latlons = coord[grep(degmin_pattern_1, coord)]
    if (length(latlons) > 0)
        {
            # Check if the pattern is present in the datafile.
            ll_lat = strsplit(latlons, "[Dd]eg|[Dd]egrees")
            deg_lat = as.numeric(sapply(ll_lat, "[[", 1))
            ll_lon = sapply(ll_lat, "[[", 2)
            # Eliminate all spaces not preceding longitudinal degrees.
            ll_lon = gsub("\\s$", "", ll_lon)
            ll_lon = gsub("^\\s", "", ll_lon)
            ll_lon = gsub("\\s[A-Za-z,;']", "", ll_lon)
            ll_lon = strsplit(ll_lon, " ")
            deg_lon = as.numeric(sapply(ll_lon, "[[", 2))
            ll_min_lat = sapply(ll_lon, "[[", 1)
            ll_min_lat = strsplit(ll_min_lat, "'|min|;| |,")
            min_lat_as_degrees = as.numeric(sapply(ll_min_lat, "[[", 1))/60
            ll_min_lon = sapply(ll_lat, "[[", 3)
            # Eliminate all spaces not preceding longitudinal degrees.
            ll_min_lon = gsub("\\s$", "", ll_min_lon)
            ll_min_lon = gsub("^\\s", "", ll_min_lon)
            ll_min_lon = gsub("\\s[A-Za-z,;']", "", ll_min_lon)
            ll_min_lon = strsplit(ll_min_lon, "'|min")
            min_lon_as_degrees = as.numeric(sapply(ll_min_lon, "[[", 1))/60
            lat = deg_lat + min_lat_as_degrees
            lon = deg_lon + min_lon_as_degrees
            lat_sign = rep("+", length(latlons))
            lon_sign = rep("+", length(latlons))
            lat_sign[grep("S", latlons)] = "-"
            lon_sign[grep("W", latlons)] = "-"
            lat[lat_sign == "-"] = -lat[lat_sign == "-"]
            lon[lon_sign == "-"] = -lon[lon_sign == "-"]
            all_lats[grep(degmin_pattern_1, coord)] = lat
            all_lons[grep(degmin_pattern_1, coord)] = lon
        }  # End if(length(latlons)>0).

    # Remove the old field 'Lat_lon' and add two different fields, one Lat (for
    # Latitude) and the Long (for Longitude), in the original dataset.
    input.new = cbind(input[, c(1:21)], Lat = all_lats, Long = all_lons, input[,
        c(23:25)])
    utils::write.table(input.new, file = output, sep = "\t", row.names = FALSE)
    return(input.new)
}  # End of the function.
