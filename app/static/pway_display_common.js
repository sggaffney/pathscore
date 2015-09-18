/**
 * Created by sgg on 2015-09-18.
 */
function create_plot_div(pway_obj, n_pathways, projId){
	var ngenes = pway_obj.geneSet.length;
	var max_height = 25*ngenes + 80;
	var pval = pway_obj['pval'];
	var pq_str = null;
	if(pval == '0.000e+00')
		pq_str = 'p<1e-16, q<' + Number(1e-16 * n_pathways).toPrecision(3);
	else
		pq_str = 'p=' + Number(pval).toPrecision(3) + ', q=' + Number(pval * n_pathways).toPrecision(3);
	div_str ='<div class="panel panel-primary">' +
		'<div class="panel-heading">' +
		'<h2 class="panel-title"><a href="'+ root_url + pway_obj['url'] +'" target="_blank">'+ pway_obj['name'] +
		'</a> (' + pq_str + ') </h2>' + pway_obj['brief'] + ' (' + pway_obj['contrib'] + ')</div>' +
		'<div class="panel-body" style="padding:0"><code class="p-lengths">Gene lengths: ' +
		'<span class="keyword-text">min</span>=' + pway_obj['lengths'][0] + 'kbp (' + pway_obj['lengths'][1] + '); ' +
		'<span class="keyword-text">max</span>=' + pway_obj['lengths'][2] + 'kbp (' + pway_obj['lengths'][3] + '); ' +
		'<span class="keyword-text">avg</span>=' + pway_obj['lengths'][4] + 'kbp; ' +
		'<span class="keyword-text">var</span>=' + pway_obj['lengths'][5] + driver_text + '</code>'
		+'<div class="svg_target ui-widget-content">' +
		'<img src="static/data/{{ user_id }}/' + projId + '/pathways_svg/' + pway_obj['id'] + '.svgz" ' +
		'alt="' + pway_obj['id'] + '.svgz" ></img>'
		+'</div>'
		+'<div class="svg_matrix ui-widget-content"><div class="pway_size" style="display:none">'
		+ JSON.stringify({"n_genes":pway_obj.geneSet.length, "expanded":false}) + '</div>' +
		'<img src="static/data/{{ user_id }}/' + projId + '/matrix_svg/' + pway_obj['id'] + '.svgz" ' +
		'px alt="' + pway_obj['id'] + '.svgz" style="max-height:' + max_height + 'px;"></img>' //
		+'</div></div>';
	return div_str;
}
