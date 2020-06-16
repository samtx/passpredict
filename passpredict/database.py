# create SQLAlchemy declarative base

import sqlalchemy
import sqlalchemy.ext.declarative

SQLALCHEMY_DATABASE_URL = "sqlite:///./sqlite.db"

engine = sqlalchemy.create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args={'check_same_thread': False}
)

Base = sqlalchemy.ext.declarative.declarative_base()


