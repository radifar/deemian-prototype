from pydantic import BaseModel, Field, validator


class Interaction(BaseModel):
    """Contain the specification for interaction identification"""

    measure: str


class Presentation(BaseModel):
    """Contain the specification for presentation format and granularity"""

    granularity: str
    field: list
    output_file: str = Field(alias="output file")

    @validator("field", pre=True)
    def split_field(cls, value):
        if isinstance(value, str):
            field_list = value.split(",")
            field_list = [field.strip() for field in field_list]
            return field_list
        elif isinstance(value, list):
            return value


class Specification(BaseModel):
    """Contain the root specification"""

    load_molecule: dict = Field(alias="load molecule")
    selection: dict
    interaction: Interaction
    presentation: Presentation
