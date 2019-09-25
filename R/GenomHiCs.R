#' GenomHiCs class.
#'
#' @slot query_regions query_regions.
#' @slot contacts Contacts.
#' @slot annotations Annotations
#' @slot structures Structures
#'
#' @name GenomHiCs-class
#' @rdname GenomHiCs-class
setClass("GenomHiCs",
         slots=list(query_regions="list",
                    contacts="GInteractions", 
                    annotations="GRangesList",
                    structures="GRangesList"))

#' Add an annotation to a GenomHiCs object.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
query_regions <- function(x) {
    return(x@query_regions)
}

#' Add an annotation to a GenomHiCs object.
#' @param x A GenomHiCs object.
#' @export
contacts <- function(x) {
    return(x@contacts)
}

#' Add an annotation to a GenomHiCs object.
#' @param x A GenomHiCs object.
#' @export
annotations <- function(x) {
    return(x@annotations)
}

#' Add an annotation to a GenomHiCs object.
#' @param x A GenomHiCs object.
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

#' Build a GenomHiCs object.
#' @export
GenomHiCs <- function(query_regions, contacts, annotations, structures=NULL) {
     res = methods::new("GenomHiCs",
                        query_regions=query_regions,
                        contacts=contacts,
                        annotations=annotations,
                        structures=structures)

    for(annot_name in names(annotations)) {
        res = add_annotation(res, annotations[[annot_name]], annot_name)
    }

    res
}

#' Add an annotation to a GenomHiCs object.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
add_annotation <- function(x, annotation, annot_name) {
    regions(x@contacts) = annotate_count(regions(x@contacts), 
                                         annotation, 
                                         annot_name)
                                       
    for(q_name in names(x@query_regions)) {
        x@query_regions[[q_name]] = annotate_count(x@query_regions[[q_name]], 
                                                   annotation, 
                                                   annot_name)
    }
    
    mcols(x@structures, level="within")[,annot_name] = 0
    for(s_name in names(x@structures)) {
        x@structures[[s_name]] = annotate_count(x@structures[[s_name]],
                                                annotation, 
                                                annot_name)
    }
    
    x
}

#' Determine if query_regions are window-bound by the given annotation.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_window_bound <- function(x, query_name, annot_name) {
    # 3a. Identify genes with binding in their window.
    gene_overlaps = InteractionSet::findOverlaps(x@query_regions[[query_name]], regions(x@contacts))

    win_bound = rep(FALSE, length(x@query_regions[[query_name]]))
    hit_binding = mcols(regions(x@contacts)[subjectHits(gene_overlaps)])[[annot_name]] > 0
    hit_df = data.frame(query=queryHits(gene_overlaps),
                        subject=subjectHits(gene_overlaps),
                        binding=hit_binding)
                   
    bind_status = hit_df %>% 
                    dplyr::group_by(query) %>% 
                    dplyr::summarize(Bound=any(binding))
    win_bound[bind_status$query] = bind_status$Bound
    
    win_bound
}

#' Determine if query_regions are contact-bound by the given annotation.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_contact_bound <- function(x, query_name, annot_name) {
    # 3b. Identify genes with GR binding in contacts.
    promoter_hic_overlap = InteractionSet::findOverlaps(x@contacts, x@query_regions[[query_name]])
    
    f_a = anchors(x@contacts, type="first")
    s_a = anchors(x@contacts, type="second")
    has_interaction = mcols(f_a)[[annot_name]] | mcols(s_a)[[annot_name]]
    
    contact_bound = rep(FALSE, length(x@query_regions[[query_name]]))
    contact_bound[S4Vectors::subjectHits(promoter_hic_overlap)] = has_interaction[S4Vectors::queryHits(promoter_hic_overlap)]
    
    contact_bound
}

#' Determine if query_regions are structure-bound by the given annotation.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
is_structure_bound = function(x, query_name, annot_name) {
    res = list()
    for(s_name in names(x@structures)) {
        s_overlap = GenomicRanges::findOverlaps(x@query_regions[[query_name]], x@structures[[s_name]])
        res[[s_name]] = rep(FALSE, length(x@query_regions[[query_name]]))
        res[[s_name]][queryHits(s_overlap)] = mcols(x@structures[[s_name]])[[annot_name]][subjectHits(s_overlap)] > 0
    }

    res
}

#' Retrieves all binding types for the given annotation.
#' @param x A GenomHiCs object.
#' @export
distant_binding = function(x, query_name, annot_name) {
    results=list(Window=is_window_bound(x, query_name, annot_name),
                 Contact=is_contact_bound(x, query_name, annot_name),
                 Structure=is_structure_bound(x, query_name, annot_name))
                 
    results
}

#' Retrieves all binding types for the given annotation and annotates the query_regions.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors mcols
#' @export
annotate_distant_binding <- function(x, query_name, annot_name) {
    bind_info = distant_binding(x, query_name, annot_name)
    mcols(x@query_regions[[query_name]])[[paste0(annot_name, "_", "Window")]] = bind_info$Window
    mcols(x@query_regions[[query_name]])[[paste0(annot_name, "_", "Contact")]] = bind_info$Contact
    
    for(s_name in names(bind_info$Structure)) {
        s_col = paste0(annot_name, "_", s_name)
        mcols(x@query_regions[[query_name]])[[s_col]] = bind_info$Structure[[s_name]]
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

#' Get a list of site directly binding the query_regions.
#' @param x A GenomHiCs object.
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @export
get_direct_sites <- function(x, query_name, annot_name) {
    # Direct binding
    promoter_overlap = GenomicRanges::findOverlaps(x@query_regions[[query_name]], x@annotations[[annot_name]])
    build_site_list(x@query_regions[[query_name]],
                    S4Vectors::queryHits(promoter_overlap),
                    x@annotations[[annot_name]],
                    S4Vectors::subjectHits(promoter_overlap))
}

#' Get a list of site binding within the windows of the query_regions.
#' @param x A GenomHiCs object.
#' @importFrom dplyr inner_join
#' @export
get_window_sites <- function(x, query_name, annot_name) {
    # In-window binding
    gene_overlaps = InteractionSet::findOverlaps(x@query_regions[[query_name]], regions(x@contacts))
    window_overlap = GenomicRanges::findOverlaps(regions(x@contacts), x@annotations[[annot_name]])
    overlap_df = dplyr::inner_join(as.data.frame(gene_overlaps), 
                                   as.data.frame(window_overlap),
                                   c(subjectHits="queryHits"))
    build_site_list(x@query_regions[[query_name]],
                    overlap_df$queryHits,
                    x@annotations[[annot_name]],
                    overlap_df$subjectHits.y)
}

#' Get a list of site distantly-bound to the query_regions.
#' @param x A GenomHiCs object.
#' @importFrom InteractionSet linkOverlaps
#' @export
get_distant_sites <- function(x, query_name, annot_name) {
    # Distant binding
    promoter_distant = InteractionSet::linkOverlaps(x@contacts, 
                                                    x@query_regions[[query_name]], 
                                                    x@annotations[[annot_name]])
    build_site_list(x@query_regions[[query_name]],
                    promoter_distant$subject1,
                    x@annotations[[annot_name]],
                    promoter_distant$subject2)
}

#' Get a list of all sites connected to the query_regions.
#' @param x A GenomHiCs object.
#' @importFrom GenomicRanges GRangesList
#' @export
get_all_sites <- function(x, query_name, annot_name) {
    direct_sites = get_direct_sites(x, query_name, annot_name)
    window_sites = get_window_sites(x, query_name, annot_name)
    distant_sites = get_distant_sites(x, query_name, annot_name)

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
