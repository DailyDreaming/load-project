import uuid


def deterministic_uuid(accession_id):
    namespace_uuid = uuid.UUID('296e1002-1c99-4877-bb7f-bb6a3b752638')
    return uuid.uuid5(namespace_uuid, accession_id)
