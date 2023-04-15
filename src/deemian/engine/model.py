from pydantic import BaseModel, Field


class Interaction(BaseModel):
    """Contain the specification for interaction identification"""

    measure: str


class Presentation(BaseModel):
    """Contain the specification for presentation format and granularity"""

    granularity: str
    field: list
    output_file: str = Field(alias="output file")


class Specification(BaseModel):
    """Contain the root specification"""

    load_molecule: dict = Field(alias="load molecule")
    selection: dict
    interaction: Interaction
    presentation: Presentation
