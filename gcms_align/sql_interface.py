'''
Created on Oct. 7, 2021

@author: Ryan
'''
import pathlib
import getpass

import sqlalchemy.orm
from sqlalchemy import Column, Integer, Float, String, ForeignKey

from .sql_base import Base
from . import identification

DEBUG_SQL = False


class SQL_Interface(object):
    ''' Maintain an sqlalchemy session and handle any objects that need conversion to upload '''

    def __init__(self):
        self.engine = None
        self.session = None

    def link_sqlite(self, file_path: pathlib.Path, clear_file: bool =True):
        self.engine = sqlalchemy.create_engine(r"sqlite:///"+str(file_path), 
                                               connect_args={"check_same_thread": False},
                                               echo=DEBUG_SQL, future=True)
        self.session = sqlalchemy.orm.Session(self.engine)
        if clear_file:
            Base.metadata.drop_all(self.engine)
        Base.metadata.create_all(self.engine)  # Create empty tables in the file

    def link_cloud_postgre(self, project_input, clear_database: bool=False):
        ''' Remote db setup '''
        input_cfg = project_input.config
        self.engine = sqlalchemy.create_engine(f"postgresql://{input_cfg.db_username}:{input_cfg.db_password}@{input_cfg.remote_db_location}:5432/{input_cfg.db_name}",
                                               connect_args={"check_same_thread": False})
        self.session = sqlalchemy.orm.Session(self.engine)
        if clear_database:
            Base.metadata.drop_all(self.engine)
        Base.metadata.create_all(self.engine)  # Create empty tables in the file

    def update_project_table(self, project_input, analysis_set):
        source_settings_storage = project_input.config
        project_input.settings_sql = AnalysisSettings_SQL(project=project_input, version=source_settings_storage.settings_version)
        project_input.settings_sql.save_attributes(self._convert_settings(source_settings_storage))
        for aligned_entry in analysis_set.alignment_data_groups:
            aps = aligned_entry.sql_ref
            project_input.aligned_peaks.append(aps)
        self.session.add(project_input)
        self.commit()

    def _convert_settings(self, input_object) -> dict:
        key_list = [x for x in dir(input_object) if not x.startswith("__")]
        result = {}
        for key in key_list:
            result[key] = getattr(input_object, key)
        result.update(vars(input_object))
        return result

    def create_search_compounds(self):  # , target_formula: str, target_cas: str) -> typing.Optional['Compound_SQL']:
        ''' Create an object to be searched for compound information '''
        result = self.session.execute(sqlalchemy.select(identification.UserEntry))  # .where((Compound_SQL.formula == target_formula) & (Compound_SQL.cas == target_cas)))
        return result.fetchall()

    def commit(self):
        self.session.commit()


class AnalysisSettings_SQL(Base):
    __tablename__ = "Analysis_Settings"
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("Minerva_Project.id"))
    project = sqlalchemy.orm.relationship("Project", back_populates="analysis_settings")
    change_date = Column(sqlalchemy.DateTime, nullable=True)
    version = Column(String, nullable=True)

    settings_list = sqlalchemy.orm.relationship("AnalysisAttribute_SQL", back_populates="settings_parent")

    def save_attributes(self, attributes: dict):
        for key, value in attributes.items():
            if callable(value) or key.startswith("db_") or key == "_store":
                continue
            try:
                input_type = type(value).__name__
            except Exception as _e:
                input_type = "unknown"
            new_attribute = AnalysisAttribute_SQL(settings_parent=self, name=key, value=str(value), type=str(input_type))
            self.settings_list.append(new_attribute)
        user_attribute = AnalysisAttribute_SQL(settings_parent=self, name="user_account_name", value=str(getpass.getuser()), type="str")
        self.settings_list.append(user_attribute)


class AnalysisAttribute_SQL(Base):
    __tablename__ = "Analysis_Attribute"
    id = Column(Integer, primary_key=True)
    settings_id = Column(Integer, ForeignKey("Analysis_Settings.id"))
    settings_parent = sqlalchemy.orm.relationship("AnalysisSettings_SQL", back_populates="settings_list")
    name = Column(String, nullable=False)
    value = Column(String, nullable=True)
    type = Column(String, nullable=False)


class AlignedPeakSet_SQL(Base):
    ''' Required due to inheritance conflicts '''
    __tablename__ = "Aligned_Peak_Set"
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("Minerva_Project.id"))
    project = sqlalchemy.orm.relationship("Project", back_populates="aligned_peaks")
    gc_peaks = sqlalchemy.orm.relationship("Peak", back_populates="aligned_set")
    compound_match = sqlalchemy.orm.relationship("CompoundMatch_SQL", back_populates="aligned_set")

    average_rt = Column(Float, nullable=False)
    total_area = Column(Float, nullable=False)
    count = Column(Integer, nullable=False)


class CompoundMatch_SQL(Base):
    ''' Map table: Links compounds to aligned sets, not directly Peak '''
    __tablename__ = "Compound_Match_Map"
    id = Column(Integer, primary_key=True)
    aligned_id = Column(Integer, ForeignKey("Aligned_Peak_Set.id"))
    aligned_set = sqlalchemy.orm.relationship("AlignedPeakSet_SQL", back_populates="compound_match")
    compound_id = Column(Integer, ForeignKey("User_ID_Compound.id"))
    compound = sqlalchemy.orm.relationship("UserEntry", back_populates="aligned_match")

    score = Column(Float, nullable=True)
    source = Column(String, nullable=True)
    match_type = Column(String, nullable=True)
