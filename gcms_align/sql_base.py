'''
Created on Oct. 28, 2022

@author: Ryan
'''
import sqlalchemy.orm


mapper_registry = sqlalchemy.orm.registry()
Base = mapper_registry.generate_base()  # alternatively: orm.declarative_base

