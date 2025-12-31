"""Repository for loading Gene domain objects."""

from typing import Dict, List, Optional, Set

from sqlalchemy import text

from ..domain import Gene, GeneSet, SizeAlgorithm
from .. import db


class GeneRepository:
    """
    Repository for loading gene data from the database.

    This class provides methods to load genes with their properties
    and calculate background genome sizes.

    Attributes:
        _gene_cache: Cache of loaded Gene objects by entrez_id
        _bmr_table: Optional custom BMR table name
    """

    def __init__(self, bmr_table: Optional[str] = None):
        """
        Initialize the repository.

        Args:
            bmr_table: Optional custom BMR table name for mutation rates.
                      If None, uses the default refs.entrez_length table.
        """
        self._gene_cache: Dict[int, Gene] = {}
        self._bmr_table = bmr_table
        self._table_name = bmr_table if bmr_table else 'refs.entrez_length'

    def get_gene(self, entrez_id: int) -> Optional[Gene]:
        """
        Get a single gene by entrez ID.

        Args:
            entrez_id: NCBI Entrez gene ID

        Returns:
            Gene object or None if not found
        """
        if entrez_id in self._gene_cache:
            return self._gene_cache[entrez_id]

        genes = self.get_genes({entrez_id})
        return genes.get(entrez_id)

    def get_genes(self, entrez_ids: Set[int]) -> Dict[int, Gene]:
        """
        Load multiple genes by entrez ID.

        Args:
            entrez_ids: Set of entrez IDs to load

        Returns:
            Dict mapping entrez_id to Gene object
        """
        # Check cache first
        result = {}
        missing = set()
        for eid in entrez_ids:
            if eid in self._gene_cache:
                result[eid] = self._gene_cache[eid]
            else:
                missing.add(eid)

        if not missing:
            return result

        # Load from database
        ids_str = ','.join(str(i) for i in missing)
        cmd = text(f"""
            SELECT l.entrez_id, l.hugo_symbol, l.length_bp, l.per_Mb, l.effective_bp
            FROM {self._table_name} l
            WHERE l.entrez_id IN ({ids_str})
        """)

        rows = db.session.execute(cmd).all()
        for row in rows:
            gene = Gene(
                entrez_id=int(row[0]),
                symbol=row[1],
                length_bp=int(row[2]) if row[2] else 0,
                mutation_rate=float(row[3]) if row[3] else 1.0
            )
            self._gene_cache[gene.entrez_id] = gene
            result[gene.entrez_id] = gene

        return result

    def get_all_genes(self) -> GeneSet:
        """
        Load all genes in the reference database.

        Returns:
            GeneSet containing all genes
        """
        cmd = text(f"""
            SELECT entrez_id, hugo_symbol, length_bp, per_Mb, effective_bp
            FROM {self._table_name}
        """)

        genes = []
        rows = db.session.execute(cmd).all()
        for row in rows:
            gene = Gene(
                entrez_id=int(row[0]),
                symbol=row[1],
                length_bp=int(row[2]) if row[2] else 0,
                mutation_rate=float(row[3]) if row[3] else 1.0
            )
            self._gene_cache[gene.entrez_id] = gene
            genes.append(gene)

        return GeneSet(genes)

    def get_background_size(
        self,
        algorithm: SizeAlgorithm,
        exclude_symbols: Optional[Set[str]] = None
    ) -> int:
        """
        Get background genome size for the given algorithm.

        Args:
            algorithm: Size calculation algorithm
            exclude_symbols: Optional set of gene symbols to exclude

        Returns:
            Background size (gene count, bp sum, or effective bp sum)
        """
        if algorithm == SizeAlgorithm.GENE_COUNT:
            field = "count(*)"
        elif algorithm == SizeAlgorithm.GENE_LENGTH:
            field = "sum(length_bp)"
        elif algorithm == SizeAlgorithm.BMR_LENGTH:
            field = "sum(effective_bp)"
        else:
            raise ValueError(f"Unknown algorithm: {algorithm}")

        exclude_clause = ""
        if exclude_symbols:
            symbols_str = ','.join(f"'{s}'" for s in exclude_symbols)
            exclude_clause = f"WHERE hugo_symbol NOT IN ({symbols_str})"

        cmd = text(f"SELECT {field} FROM {self._table_name} {exclude_clause}")
        result = db.session.execute(cmd).scalar()
        return int(result) if result else 0

    def get_genes_in_pathways(self) -> GeneSet:
        """
        Load only genes that are members of at least one pathway.

        Returns:
            GeneSet containing genes with pathway membership
        """
        cmd = text(f"""
            SELECT DISTINCT l.entrez_id, l.hugo_symbol, l.length_bp, l.per_Mb
            FROM {self._table_name} l
            INNER JOIN refs.pathway_gene_link pgl ON l.entrez_id = pgl.entrez_id
        """)

        genes = []
        rows = db.session.execute(cmd).all()
        for row in rows:
            gene = Gene(
                entrez_id=int(row[0]),
                symbol=row[1],
                length_bp=int(row[2]) if row[2] else 0,
                mutation_rate=float(row[3]) if row[3] else 1.0
            )
            self._gene_cache[gene.entrez_id] = gene
            genes.append(gene)

        return GeneSet(genes)

    def clear_cache(self):
        """Clear the gene cache."""
        self._gene_cache.clear()
