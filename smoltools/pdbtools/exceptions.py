class ChainNotFound(KeyError):
    def __init__(self, structure_id: str, model_id: str, chain_id: str):
        message = f'Chain {structure_id}/{model_id}/{chain_id} not in structure'
        super().__init__(message)

class NoResiduesFound(ValueError):
    def __init__(self) -> None:
        message = 'No residues matching filter criteria found.'
        super().__init__(message)


class NoAtomsFound(ValueError):
    def __init__(self) -> None:
        message = 'No atoms matching filter criteria found.'
        super().__init__(message)
