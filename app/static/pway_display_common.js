/**
 * Created by sgg on 2015-09-18.
 */
function create_plot_div(pway_obj, n_pathways, user_id, projId, matrix_dir, n_genes, opt_header){
    matrix_dir = matrix_dir || "matrix_svg";
    n_genes = n_genes || pway_obj.geneSet.length;
    opt_header = opt_header || "";
    var driver_text = get_driver_str(pway_obj);
	var max_height = 25 * n_genes + 80;
	var pval = pway_obj['pval'];
	var pq_str = null;
	if(pval == '0.000e+00')
		pq_str = '<em>P</em><1e-16, <em>Q</em><' + Number(1e-16 * n_pathways).toPrecision(1);
	else {
        pq_str = '<em>P</em>=' + get_nice_precision(pval) + ', <em>Q</em>=' + get_nice_precision(pval * n_pathways);
    }
	div_str ='<div class="panel panel-primary">' +
		'<div class="panel-heading">' + opt_header +
		'<h2 class="panel-title"><a href="'+ root_url + pway_obj['url'] +'" target="_blank">'+ pway_obj['name'] +
		'</a> (' + pq_str + ') </h2>' + pway_obj['brief'] + ' (' + pway_obj['contrib'] + ')</div>' +
		'<div class="panel-body" style="padding:0"><code class="p-lengths">Gene lengths: ' +
		'<span class="keyword-text">min</span>=' + pway_obj['lengths'][0] + 'kbp (' + pway_obj['lengths'][1] + '); ' +
		'<span class="keyword-text">max</span>=' + pway_obj['lengths'][2] + 'kbp (' + pway_obj['lengths'][3] + '); ' +
		'<span class="keyword-text">avg</span>=' + pway_obj['lengths'][4] + 'kbp; ' +
		'<span class="keyword-text">var</span>=' + pway_obj['lengths'][5] + driver_text + '</code>'
		+'<div class="svg_target ui-widget-content">' +
		'<img src="static/data/' + user_id + '/' + projId + '/pathways_svg/' + pway_obj['id'] + '.svgz" ' +
		'alt="' + pway_obj['id'] + '.svgz" ></img>'
		+'</div>'
		+'<div class="svg_matrix ui-widget-content"><div class="pway_size" style="display:none">'
		+ JSON.stringify({"n_genes":pway_obj.geneSet.length, "expanded":false}) + '</div>' +
		'<img src="static/data/' + user_id + '/' + projId + '/' + matrix_dir + '/' + pway_obj['id'] + '.svgz" ' +
		'px alt="' + pway_obj['id'] + '.svgz" style="max-height:' + max_height + 'px;"></img>' //
		+'</div></div>';
	return div_str;
}

function get_driver_str(pway_obj){
    vogelstein = ['ABL1', 'ACVR1B', 'AKT1', 'ALK', 'APC', 'AR', 'ARID1A', 'ARID1B', 'ARID2', 'ASXL1', 'ATM', 'ATRX', 'AXIN1', 'B2M', 'BAP1', 'BCL2', 'BCOR', 'BRAF', 'BRCA1', 'BRCA2', 'CARD11', 'CASP8', 'CBL', 'CDC73', 'CDH1', 'CDKN2A', 'CEBPA', 'CIC', 'CREBBP', 'CRLF2', 'CSF1R', 'CTNNB1', 'CYLD', 'DAXX', 'DNMT1', 'DNMT3A', 'EGFR', 'EP300', 'ERBB2', 'EZH2', 'AMER1', 'FBXW7', 'FGFR2', 'FGFR3', 'FLT3', 'FOXL2', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'H3F3A', 'HIST1H3B', 'HNF1A', 'HRAS', 'IDH1', 'IDH2', 'JAK1', 'JAK2', 'JAK3', 'KDM5C', 'KDM6A', 'KIT', 'KLF4', 'KRAS', 'MAP2K1', 'MAP3K1', 'MED12', 'MEN1', 'MET', 'MLH1', 'KMT2B', 'KMT2C', 'MPL', 'MSH2', 'MSH6', 'MYD88', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NOTCH1', 'NOTCH2', 'NPM1', 'NRAS', 'PAX5', 'PBRM1', 'PDGFRA', 'PHF6', 'PIK3CA', 'PIK3R1', 'PPP2R1A', 'PRDM1', 'PTCH1', 'PTEN', 'PTPN11', 'RB1', 'RET', 'RNF43', 'RUNX1', 'SETD2', 'SETBP1', 'SF3B1', 'SMAD2', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMO', 'SOCS1', 'SOX9', 'SPOP', 'SRSF2', 'STAG2', 'STK11', 'TET2', 'TNFAIP3', 'TRAF7', 'TP53', 'TSC1', 'TSHR', 'U2AF1', 'VHL', 'WT1', 'ALPK2', 'BCLAF1', 'CD1D', 'CEP76', 'CHD8', 'DNER', 'ELF3', 'EZH1', 'ARHGAP35', 'HIST1H4E', 'HLA-B', 'IRF6', 'MAP4K3', 'MBD1', 'MGA', 'MYOCD', 'PCBP1', 'QKI', 'RAD21', 'RHEB', 'RHOA', 'RPL5', 'RXRA', 'SETDB1', 'SOS1', 'STX2', 'TAP1', 'TNF', 'TP53BP1', 'TPX2', 'TRIM23', 'ZNF750', 'ZRANB3'];
    var driver_text = "";
    var driver_list = [];
    $.each(pway_obj.geneSet, function(gi, gene_name){
        if( $.inArray(gene_name, vogelstein)>-1 ){
            driver_list.push(gene_name)
        }
    });
    driver_list.sort();
    if (driver_list.length){
        driver_text = '<br>&nbsp;Vogelstein drivers: ' + driver_list.join(', ');
    }
    return driver_text;
}

function get_nice_precision(value, precision){
    precision = precision || 2;
    var presVal = Number(value).toPrecision(precision);
    if(value >= 1)
        return 1;
    else if(value >= 0.001)
        return presVal;
    else
        return Number(presVal).toExponential();
}