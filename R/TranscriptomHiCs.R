#' TranscriptomHiCs class.
#'
#' @slot promoters Promoters.
#' @slot contacts Contacts.
#' @slot
#' @slot 
#'
#' @name TranscriptomHiCs-class
#' @rdname TranscriptomHiCs-class
setClass("TranscriptomHiCs",
         slots=list(promoters="GRanges",
                    contacts="GInteractions", 
                    annotations="GRangesList",
                    structures="GRangesList"))

#' Add an annotation to a TranscriptomHiCs object.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
promoters <- function(x) {
    return(x@promoters)
}

#' Add an annotation to a TranscriptomHiCs object.
#' @param x A TranscriptomHiCs object.
#' @export
contacts <- function(x) {
    return(x@contacts)
}

#' Add an annotation to a TranscriptomHiCs object.
#' @param x A TranscriptomHiCs object.
#' @export
annotations <- function(x) {
    return(x@annotations)
}

#' Add an annotation to a TranscriptomHiCs object.
#' @param x A TranscriptomHiCs object.
#' @export
structures <- function(x) {
    return(x@structures)
}

#' Adds an annotation count column to a GRanges object.
#' @importFrom S4Vectors mcols
#' @keywords internal
annotate_count = function(query_regions, annotation, label) {
    overlap_res = GenomicRanges::countOverlaps(query_regions, annotation)
    mcols(query_regions)[[label]] = overlap_res    

    query_regions
}

#' Build a TranscriptomHiCs object.
#' @export
transcriptomHiCs <- function(promoters, contacts, annotations, structures=NULL) {
     res = methods::new("TranscriptomHiCs",
                        promoters=promoters,
                        contacts=contacts,
                        annotations=annotations,
                        structures=structures)

    for(annot_name in names(annotations)) {
        res = add_annotation(res, annotations[[annot_name]], annot_name)
    }

    res
}

#' Add an annotation to a TranscriptomHiCs object.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
add_annotation <- function(x, annotation, annot_name) {
    regions(x@contacts) = annotate_count(regions(x@contacts), 
                                         annotation, 
                                         annot_name)
                                       
    x@promoters = annotate_count(x@promoters, 
                                 annotation, 
                                 annot_name)
    
    mcols(x@structures, level="within")[,annot_name] = 0
    for(s_name in names(x@structures)) {
        x@structures[[s_name]] = annotate_count(x@structures[[s_name]],
                                                annotation, 
                                                annot_name)
    }
    
    x
}

#' Determine if promoters are window-bound by the given annotation.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_window_bound <- function(x, colname) {
    # 3a. Identify genes with binding in their window.
    gene_overlaps = GenomicRanges::findOverlaps(x@promoters, regions(x@contacts))

    win_bound = rep(FALSE, length(x@promoters))
    hit_binding = mcols(regions(x@contacts)[subjectHits(gene_overlaps)])[[colname]] > 0
    hit_df = data.frame(query=queryHits(gene_overlaps),
                        subject=subjectHits(gene_overlaps),
                        binding=hit_binding)
                   
    bind_status = hit_df %>% 
                    dplyr::group_by(query) %>% 
                    dplyr::summarize(Bound=any(binding))
    win_bound[bind_status$query] = bind_status$Bound
    
    win_bound
}

#' Determine if promoters are contact-bound by the given annotation.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_contact_bound <- function(x, colname) {
    # 3b. Identify genes with GR binding in contacts.
    promoter_hic_overlap = GenomicRanges::findOverlaps(x@contacts, x@promoters)
    
    f_a = anchors(x@contacts, type="first")
    s_a = anchors(x@contacts, type="second")
    has_interaction = mcols(f_a)[[colname]] | mcols(s_a)[[colname]]
    
    contact_bound = rep(FALSE, length(x@promoters))
    contact_bound[S4Vectors::subjectHits(promoter_hic_overlap)] = has_interaction[S4Vectors::queryHits(promoter_hic_overlap)]
    
    contact_bound
}

#' Determine if promoters are structure-bound by the given annotation.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_structure_bound = function(x, colname) {
    res = list()
    for(s_name in names(x@structures)) {
        s_overlap = GenomicRanges::findOverlaps(x@promoters, x@structures[[s_name]])
        res[[s_name]] = rep(FALSE, length(x@promoters))
        res[[s_name]][queryHits(s_overlap)] = mcols(x@structures[[s_name]])[[colname]][subjectHits(s_overlap)] > 0
    }

    res
}

#' Retrieves all binding types for the given annotation.
#' @param x A TranscriptomHiCs object.
#' @export
distant_binding = function(x, annot_name) {
    results=list(Window=is_window_bound(x, annot_name),
                 Contact=is_contact_bound(x, annot_name),
                 Structure=is_structure_bound(x, annot_name))
                 
    results
}

#' Retrieves all binding types for the given annotation and annotates the promoters.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
annotate_distant_binding <- function(x, annot_name) {
    bind_info = distant_binding(x, annot_name)
    mcols(x@promoters)[[paste0(annot_name, "_", "Window")]] = bind_info$Window
    mcols(x@promoters)[[paste0(annot_name, "_", "Contact")]] = bind_info$Contact
    
    for(s_name in names(bind_info$Structure)) {
        s_col = paste0(annot_name, "_", s_name)
        mcols(x@promoters)[[s_col]] = bind_info$Structure[[s_name]]
    }
    
    x
}

#' Go from promoter/enhancer indices to a GRangesList of sites.
#' 
#' @importFrom GenomicRanges GRangesList
#' @keywords internal
build_site_list = function(promoter_regions, 
                           promoter_ids, 
                           query_regions,
                           query_ids) {
    unique_query = unique(promoter_ids)
    query_by_gene = lapply(unique_query, function(x) {
        query_regions[query_ids[promoter_ids==x] ]
    })
    names(query_by_gene) = names(promoter_regions)[unique_query]
    query_by_gene = GenomicRanges::GRangesList(query_by_gene)

    query_by_gene
}

#' Get a list of site directly binding the promoters.
#' @param x A TranscriptomHiCs object.
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
get_direct_sites <- function(x, query_regions) {
    # Direct binding
    promoter_overlap = GenomicRanges::findOverlaps(x@promoters, query_regions)
    build_site_list(x@promoters,
                    S4Vectors::queryHits(promoter_overlap),
                    query_regions,
                    S4Vectors::subjectHits(promoter_overlap))
}

#' Get a list of site binding within the windows of the promoters.
#' @param x A TranscriptomHiCs object.
#' @importFrom dplyr inner_join
#' @export
get_window_sites <- function(x, query_regions) {
    # In-window binding
    gene_overlaps = InteractionSet::findOverlaps(x@promoters, regions(x@contacts))
    window_overlap = GenomicRanges::findOverlaps(regions(x@contacts), query_regions)
    overlap_df = dplyr::inner_join(as.data.frame(gene_overlaps), 
                                   as.data.frame(window_overlap),
                                   c(subjectHits="queryHits"))
    build_site_list(x@promoters,
                    overlap_df$queryHits,
                    query_regions,
                    overlap_df$subjectHits.y)
}

#' Get a list of site distantly-bound to the promoters.
#' @param x A TranscriptomHiCs object.
#' @importFrom InteractionSet linkOverlaps
#' @export
get_distant_sites <- function(x, query_regions) {
    # Distant binding
    promoter_distant = InteractionSet::linkOverlaps(x@contacts, 
                                                    x@promoters, 
                                                    query_regions)
    build_site_list(x@promoters,
                    promoter_distant$subject1,
                    query_regions,
                    promoter_distant$subject2)
}

#' Get a list of all sites connected to the promoters.
#' @param x A TranscriptomHiCs object.
#' @importFrom GenomicRanges GRangesList
#' @export
get_all_sites <- function(x, query_regions) {
    direct_sites = get_direct_sites(x, query_regions)
    window_sites = get_window_sites(x, query_regions)
    distant_sites = get_distant_sites(x, query_regions)

    unique_genes = unique(c(names(direct_sites), 
                            names(window_sites), 
                            names(distant_sites)))
    
    # All binding
    site_by_gene_any = list()
    for(gene in unique_genes) {
        grl = unlist(GenomicRanges::GRangesList(c(direct_sites[[gene]],
                                                window_sites[[gene]],
                                                distant_sites[[gene]])))
        site_by_gene_any[[gene]] = grl
    }

    return(site_by_gene_any)
}
