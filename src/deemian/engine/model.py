from pydantic import BaseModel, Field, validator


class Interaction(BaseModel):
    """Contain the detail for interaction identification"""

    measure: str


class Presentation(BaseModel):
    """Contain the detail for presentation format and granularity"""

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


class Workflow(BaseModel):
    """Contain the root workflow"""

    load_molecule: dict = Field(alias="load molecule")
    selection: dict
    interaction: Interaction
    presentation: Presentation
