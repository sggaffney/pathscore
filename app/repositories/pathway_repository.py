"""Repository for loading Pathway domain objects."""

from typing import Dict, List, Optional, Set

from sqlalchemy import text

from ..domain import Gene, GeneSet, Pathway, PathwayMetadata, SizeAlgorithm
from .. import db
from .gene_repository import GeneRepository


class PathwayRepository:
    """
    Repository for loading pathway data from the database.

    This class provides methods to load pathways with their genes
    and metadata.

    Attributes:
        _pathway_cache: Cache of loaded Pathway objects by path_id
        _gene_repo: GeneRepository for loading gene data
    """

    def __init__(self, gene_repository: Optional[GeneRepository] = None):
        """
        Initialize the repository.

        Args:
            gene_repository: Optional GeneRepository to use.
                           If None, creates a new one.
        """
        self._pathway_cache: Dict[int, Pathway] = {}
        self._gene_repo = gene_repository or GeneRepository()

    def get_pathway(self, path_id: int, load_genes: bool = True) -> Optional[Pathway]:
        """
        Get a single pathway by ID.

        Args:
            path_id: Pathway identifier
            load_genes: Whether to load the genes (default True)

        Returns:
            Pathway object or None if not found
        """
        if path_id in self._pathway_cache:
            return self._pathway_cache[path_id]

        pathways = self.get_pathways({path_id}, load_genes=load_genes)
        return pathways.get(path_id)

    def get_pathways(
        self,
        path_ids: Set[int],
        load_genes: bool = True
    ) -> Dict[int, Pathway]:
        """
        Load multiple pathways by ID.

        Args:
            path_ids: Set of pathway IDs to load
            load_genes: Whether to load the genes

        Returns:
            Dict mapping path_id to Pathway object
        """
        # Check cache first
        result = {}
        missing = set()
        for pid in path_ids:
            if pid in self._pathway_cache:
                result[pid] = self._pathway_cache[pid]
            else:
                missing.add(pid)

        if not missing:
            return result

        # Load pathway metadata
        ids_str = ','.join(str(i) for i in missing)
        cmd = text(f"""
            SELECT path_id, pathway_name, info_url, description_brief, contributor
            FROM refs.pathways
            WHERE path_id IN ({ids_str})
        """)

        pathway_data = {}
        rows = db.session.execute(cmd).all()
        for row in rows:
            pathway_data[int(row[0])] = {
                'name': row[1],
                'metadata': PathwayMetadata(
                    url=row[2] or '',
                    description=row[3] or '',
                    contributor=row[4] or ''
                )
            }

        # Load genes if requested
        pathway_genes: Dict[int, Set[int]] = {}
        if load_genes:
            cmd = text(f"""
                SELECT path_id, entrez_id
                FROM refs.pathway_gene_link
                WHERE path_id IN ({ids_str})
            """)
            rows = db.session.execute(cmd).all()
            for row in rows:
                pid = int(row[0])
                eid = int(row[1])
                if pid not in pathway_genes:
                    pathway_genes[pid] = set()
                pathway_genes[pid].add(eid)

            # Load gene data
            all_entrez_ids = set()
            for genes in pathway_genes.values():
                all_entrez_ids.update(genes)
            gene_map = self._gene_repo.get_genes(all_entrez_ids)

        # Build Pathway objects
        for pid, data in pathway_data.items():
            if load_genes and pid in pathway_genes:
                genes = [gene_map[eid] for eid in pathway_genes[pid]
                        if eid in gene_map]
                gene_set = GeneSet(genes)
            else:
                gene_set = GeneSet()

            pathway = Pathway(
                path_id=pid,
                name=data['name'],
                genes=gene_set,
                metadata=data['metadata']
            )
            self._pathway_cache[pid] = pathway
            result[pid] = pathway

        return result

    def get_all_pathways(self, load_genes: bool = True) -> List[Pathway]:
        """
        Load all pathways from the database.

        Args:
            load_genes: Whether to load genes for each pathway

        Returns:
            List of all Pathway objects
        """
        # Get all pathway IDs
        cmd = text("""
            SELECT DISTINCT path_id FROM refs.pathway_gene_link
        """)
        rows = db.session.execute(cmd).all()
        path_ids = {int(row[0]) for row in rows}

        pathways = self.get_pathways(path_ids, load_genes=load_genes)
        return list(pathways.values())

    def get_pathways_for_genes(self, entrez_ids: Set[int]) -> List[Pathway]:
        """
        Get all pathways that contain any of the given genes.

        Args:
            entrez_ids: Set of entrez IDs

        Returns:
            List of Pathway objects containing at least one of the genes
        """
        if not entrez_ids:
            return []

        ids_str = ','.join(str(i) for i in entrez_ids)
        cmd = text(f"""
            SELECT DISTINCT path_id
            FROM refs.pathway_gene_link
            WHERE entrez_id IN ({ids_str})
        """)
        rows = db.session.execute(cmd).all()
        path_ids = {int(row[0]) for row in rows}

        pathways = self.get_pathways(path_ids, load_genes=True)
        return list(pathways.values())

    def get_pathway_sizes(
        self,
        algorithm: SizeAlgorithm,
        exclude_symbols: Optional[Set[str]] = None
    ) -> Dict[int, float]:
        """
        Get sizes for all pathways using the given algorithm.

        This is more efficient than loading full pathways when only
        sizes are needed.

        Args:
            algorithm: Size calculation algorithm
            exclude_symbols: Optional gene symbols to exclude

        Returns:
            Dict mapping path_id to size
        """
        table_name = 'refs.entrez_length'

        if algorithm == SizeAlgorithm.GENE_COUNT:
            field = "count(DISTINCT pgl.entrez_id)"
        elif algorithm == SizeAlgorithm.GENE_LENGTH:
            field = "sum(l.length_bp)"
        elif algorithm == SizeAlgorithm.BMR_LENGTH:
            field = "sum(l.effective_bp)"
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

        exclude_clause = ""
        if exclude_symbols:
            symbols_str = ','.join(f"'{s}'" for s in exclude_symbols)
            exclude_clause = f"AND l.hugo_symbol NOT IN ({symbols_str})"

        cmd = text(f"""
            SELECT pgl.path_id, {field}
            FROM refs.pathway_gene_link pgl
            INNER JOIN {table_name} l ON pgl.entrez_id = l.entrez_id
            WHERE 1=1 {exclude_clause}
            GROUP BY pgl.path_id
        """)

        result = {}
        rows = db.session.execute(cmd).all()
        for row in rows:
            result[int(row[0])] = float(row[1]) if row[1] else 0.0

        return result

    def clear_cache(self):
        """Clear the pathway cache."""
        self._pathway_cache.clear()
