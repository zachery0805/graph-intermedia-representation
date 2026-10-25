#include "graph.h"
graph::graph(
		const char *fw_beg_file,
		const char *fw_csr_file,
		const char *bw_beg_file,
		const char *bw_csr_file)
{
	double tm=wtime();
	
	typedef int index_tt;
	typedef int vertex_tt;

	vert_count=fsize(fw_beg_file)/sizeof(index_tt) - 1;
	edge_count=fsize(fw_csr_file)/sizeof(vertex_tt);
    
    //fw
	FILE *file=fopen(fw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_beg_file<<" cannot open\n";
		exit(-1);
	}

	index_tt *tmp_beg_pos = new index_tt[vert_count+1];
	index_tt ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(fw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<fw_csr_file<<" cannot open\n";
		exit(-1);
	}

	vertex_tt *tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	fw_beg_pos = new index_t[vert_count+1];
	fw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		fw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		fw_csr[i]=(vertex_t)tmp_csr[i];
	delete[] tmp_beg_pos;
	delete[] tmp_csr;

    // bw

	file=fopen(bw_beg_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_beg_file<<" cannot open\n";
		exit(-1);
	}

	tmp_beg_pos = new index_tt[vert_count+1];
	ret=fread(tmp_beg_pos, sizeof(index_tt), vert_count+1, file);
	assert(ret==vert_count+1);
	fclose(file);

	file=fopen(bw_csr_file, "rb");
	if(file==NULL)
	{
		std::cout<<bw_csr_file<<" cannot open\n";
		exit(-1);
	}

	tmp_csr = new vertex_tt[edge_count];
	ret=fread(tmp_csr, sizeof(vertex_tt), edge_count, file);
	assert(ret==edge_count);
	fclose(file);

    //converting to uint32_t
	bw_beg_pos = new index_t[vert_count+1];
	bw_csr = new vertex_t[edge_count];
	
	for(index_t i=0;i<vert_count+1;++i)
		bw_beg_pos[i]=(index_t)tmp_beg_pos[i];

	for(index_t i=0;i<edge_count;++i)
		bw_csr[i]=(vertex_t)tmp_csr[i];

	delete[] tmp_beg_pos;
	delete[] tmp_csr;

	// Validate both CSRs: on-disk files written with 8-byte integers (or
	// otherwise mis-converted) used to load silently as garbage and make
	// every later algorithm hang or return nonsense. See CHANGELOG 2026-10-13.
	auto check = [&](const char* name, const index_t* beg, const vertex_t* csr) {
		bool ok = beg[0] == 0 && beg[vert_count] == edge_count;
		for (index_t i = 0; ok && i < vert_count; ++i) ok = beg[i] <= beg[i + 1];
		for (index_t i = 0; ok && i < edge_count; ++i) ok = csr[i] >= 0 && csr[i] < vert_count;
		if (!ok) {
			std::cerr << "graph: " << name << " CSR in " << fw_beg_file << " / " << bw_beg_file
			          << " is malformed (offsets not monotone from 0 to E, or neighbour ids outside [0,V)). "
			          << "Was it written with 8-byte integers? Refusing to load." << std::endl;
			exit(1);
		}
	};
	{
		off_t bw_v = fsize(bw_beg_file) / sizeof(index_tt) - 1, bw_e = fsize(bw_csr_file) / sizeof(vertex_tt);
		if (bw_v != vert_count || bw_e != edge_count) {
			std::cerr << "graph: forward CSR has " << vert_count << " vertices / " << edge_count
			          << " edges but backward CSR files hold " << bw_v << " / " << bw_e
			          << " entries; refusing to load " << bw_beg_file << std::endl;
			exit(1);
		}
	}
	check("forward", fw_beg_pos, fw_csr);
	check("backward", bw_beg_pos, bw_csr);

	std::cout<<"Graph load (success): "<<vert_count<<" verts, "
		<<edge_count<<" edges "<<wtime()-tm<<" second(s)\n";
}


graph* graph::fromEdges(index_t V, const std::vector<std::pair<int, int>>& edges) {
	graph* g = new graph();
	g->vert_count = V;
	g->edge_count = (index_t)edges.size();
	g->fw_beg_pos = new index_t[V + 1]();
	g->bw_beg_pos = new index_t[V + 1]();
	g->fw_csr = new vertex_t[edges.size() > 0 ? edges.size() : 1];
	g->bw_csr = new vertex_t[edges.size() > 0 ? edges.size() : 1];
	for (auto& e : edges) { g->fw_beg_pos[e.first + 1]++; g->bw_beg_pos[e.second + 1]++; }
	for (index_t v = 0; v < V; ++v) { g->fw_beg_pos[v + 1] += g->fw_beg_pos[v]; g->bw_beg_pos[v + 1] += g->bw_beg_pos[v]; }
	std::vector<index_t> fp(g->fw_beg_pos, g->fw_beg_pos + V), bp(g->bw_beg_pos, g->bw_beg_pos + V);
	for (auto& e : edges) { g->fw_csr[fp[e.first]++] = e.second; g->bw_csr[bp[e.second]++] = e.first; }
	return g;
}
