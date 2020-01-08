#!/usr/bin/env python
# this script makes assumptions about the Fus-Api, this api is subject to change.
# when re-writing this please attempt to use the dcp-cli tool, it shall allow for easier, and
# updated interface with Fusillade.

 # Project ACLS
import hca
import json
import requests
import os

def read_bearer_from_file(filename):
    with open(filename) as config:
        found = json.load(config)
        if found is None:
            raise RuntimeError('Unable to load config file')
            exit(-1)
        token_type = found["oauth2_token"]["token_type"]
        access_token = found["oauth2_token"]["access_token"]
        return (token_type,access_token)


fus_env = 'testing'
fus_url = f"https://auth.{fus_env}.data.humancellatlas.org"
fus_url_swagger = f'{fus_url}/swagger.json'
fus_url_api = f'{fus_url}/v1'
protected_projects = ["aa372dea-8469-5b80-9007-18c16a21655d"] # random project
open_projects = [
'6ae3cbfe-200e-5c03-a74a-edd266b8182b',
'06f8848d-9c54-5829-92d3-d334809ad1e2',
'a0b6322d-1da3-5481-8768-84227ad4dd1e',
'b0f40b69-943f-5959-9457-c8e53c2d480e',
'9fc2d285-804a-5989-956f-1843a0f11673',
'099c02da-23b2-5748-8618-92bc6770dc51',
'dd761426-cc9c-5fbd-8bec-93a4cc4eb999',
'096b7311-2bf7-5e61-9afb-d65c24a71243',
'7dd4e06c-c889-511c-a4d4-45b74088caa8',
'd36952f4-cfa7-5d03-b4b6-db2c31dd41c6',
'04ba7269-1301-5758-8f13-025565326f66',
'd9258dc7-985e-533a-9fc2-3ad9cc7e32ca',
'f10b6bef-febb-58dd-83ee-1180d076e53f',
'efddbf1e-a5cf-5e61-a38f-6a9f228e07c2',
'cc2112b7-9df1-5910-a7c6-6e41203130fa',
'd136eea6-03f9-5f02-86c7-c677b4c80164',
'b26137d3-a709-5492-aa74-0d783e6b628b',
'227f5c51-389c-576c-b4d3-e4da53b89f79',
'bfddbefc-f2fd-5815-89a9-a94ed667be82',
'57c9a7c8-5dc5-551b-b4af-7d37d8a87f64',
'1cafb09c-e0dc-536b-b166-5cb6debfc3cf',
'b4b128d5-61e5-510e-9d91-a151b94fbb99',
'ef6f570f-a991-5528-9649-fbf06e6eb896',
'06917f50-92aa-5e58-8376-aae1d888e8b7',
'36a7d62a-57ae-59af-ae8d-7dd2a8f1422e',
'bdfe2399-b8a6-5b6a-9f0a-a5fd81d08ff4',
'31d48835-7d9f-52ae-8cdc-ae227b63dd2c',
'53fe8e51-1646-5f68-96ac-e4d4fde67b93',
'061ec9d5-9acf-54db-9eee-555136d5ce41',
'ed008b9b-0039-5ec2-a557-f082d4ba1810',
'458cbaeb-c4b7-5537-b1f3-a5d537478112',
'389ad9f9-4a14-5a3d-b971-45dc3baf95f1',
'56483fc6-ab20-5495-bb93-8cd2ce8a322a',
'1a72144b-f4e8-5fd5-b46e-6a40eee8b6e6',
'62437ea1-3d06-5f22-b9de-a7d934138dd5',
'427157c0-993f-56ae-9e40-8cfe40ef81c5',
'11d48902-5824-5520-836b-7fc78ed02a61',
'0a8f2289-5862-5bf0-8c27-0885453de788',
'86963b4f-1e8e-5691-9ba3-465f3a789428',
'73fd591e-2310-5983-ba8a-8079c0d0b758',
'2ff83170-ec52-5b24-9962-161c558f52ba',
'd26d2ae7-4355-5ac1-8476-2e514973097e',
'd9117a4f-36e0-5912-b8cd-744a0c5306c7',
'069198f7-c2b5-5d39-988c-1cb12db4f28a',
'c2e2302f-4077-5394-9ee2-78a0ec94cbb7',
'2cdd0744-6422-57fd-8fdd-9ac2bb8bf257',
'df1875a9-1a6a-58e0-8fd2-dac0ebb3b1b2',
'7eedae3a-b350-5a3a-ab2d-8dcebc4a37b2',
'3fe16b18-e782-542b-b308-de9b26e7f69c',
'56d9146d-bc73-5327-9615-05931f1863f6',
'2b4c411d-35b1-5f97-9d6a-c9331c7f679a',
'aa372dea-8469-5b80-9007-18c16a21655d',
'f0caef8c-f839-539d-aa17-61fe04e6d3dd',
'06a318d9-54d8-5e41-aab5-f2d682fba690',
'6fb6d88c-7023-53fb-967b-ef95b2f6f5a0',
'ff0a4d85-a1c7-571c-97ae-d964eee7ecad',
'80ad934f-66ed-5c21-8f9a-7d3b0f58bcab',
'6b892786-989c-5844-b191-6d420e328fdf',
'682f2474-f875-5e0a-bf99-2b102c8c6193',
'cac7f9f2-0592-5617-9530-f63803c49f8b',
'789850ec-3540-5023-9767-fb8a4d2a21fc',
'932ae148-c3c2-5b07-91c0-2083cafe0dc1',
'ebb8c1be-6739-57b8-9ce3-aa67caa900b4',
'30fb622d-6629-527e-a681-6e2ba143af3d',
'bd4ebaac-7bcb-5069-b3ea-3a13887092e8',
'1a7ccf4f-a500-5aa6-b31a-2fff80cf8f08',
'51a21599-a014-5c5a-9760-d5bdeb80f741',
'4f5c0011-416d-5e8e-8eb6-f7cb5b0140a5',
'e8579e71-7472-5671-85b0-9841a4d06d5a',
'c366f1f5-27aa-5157-a142-110e492a3e52',
'5cd871a3-96ab-52ed-a7c9-77e91278c13d',
'39bfc05a-44ca-507a-bbf5-156bd35c5c74',
'b5a0936b-a351-54ac-8d7d-0af6926e0bdc',
'0f4d7e06-5f77-5614-8cd6-123f555dc9b1',
'b48f6f16-1b5a-5055-9e14-a8920e1bcaad',
'7acdc227-c543-5b0c-8bd8-c6fa4e30310a',
'b55c0638-d86b-5665-9ad1-0d45b937770a',
'8d1bf054-faad-5ee8-a67e-f9b8f379e6c3',
'9d65c4d0-c048-5c4f-8278-85dac99ea2ae',
'60ec348b-ff28-5d47-b0d6-b787f1885c9c',
'99ab18ff-ca15-5d7d-9c2d-5c0b537fb1c2',
'5016f45e-b86a-57ce-984e-a50605641d08',
'88564eae-cceb-5eee-957d-5f0a251fb177',
'1a85f2da-5aaa-5bc2-a6ea-9584742118e6',
'ba75d697-8712-5e75-b35b-fd3f2b66cae5']
resource_type_name = "projects"
open_resource_policy_name = "open_access"
dss_allowed_roles = "dss_admin"
token_type, access_token = read_bearer_from_file(f"{os.getenv('HOME')}/.config/hca/config.json") # ensure that you are logged in

auth_headers = {"Authorization": f'{token_type} {access_token}',
                "Content-Type": "application/json",
                "accept": "*/*"}

# Test out auth
#  If you get 401 errors, logout and login with the `hca dss logout && hca dss login`
resp = requests.get(url=f'{fus_url_api}/users',headers=auth_headers)
assert resp.status_code is requests.codes.ok

resp = requests.get(f"{fus_url_api}/resources",headers=auth_headers)
resp_data = resp.json().get("resources")
print(resp_data)
if resource_type_name in resp_data:
    print("It appears projects are already defined in Fus, exiting")
    # exit(1)

# else its not defined and we continue....
type_data = {"actions": ["dss:GetBundle", "dss:GetFile"]}
resp = requests.post(url=f"{fus_url_api}/resource/{resource_type_name}",
                     headers=auth_headers,
                     json=type_data)
print(resp)


access_policy = {
  "policy": {
      "Statement": [
           {
               "Sid": open_resource_policy_name,
               "Effect": "Allow",
               "Action": type_data['actions'],
               "Resource": [f"arn:dcp:fus:us-east-1:{fus_env}:{resource_type_name}/*"]
           }
        ]
    },
  "type": "ResourcePolicy"
}

resp = requests.post(url=f'{fus_url_api}/resource/{resource_type_name}/policy/{open_resource_policy_name}',
                     headers=auth_headers,
                     json=access_policy)
print(json.dumps(access_policy,indent=4))
print(f'{fus_url_api}/resource/{resource_type_name}/policy/{open_resource_policy_name}')
print(resp)


# for project in protected_projects:
#     resp = requests.post(url=f'{fus_url_api}/resource/{resource_type_name}/id/{project}',
#                          headers=auth_headers)
#     print(resp)
#
#     member_policy = [
#         {
#             "Add or update access":{
#                 "value":{
#                     "member_type": "dss_admin",
#                     "access_level": locked_resource_policy_name
#                 }
#              }
#         }
#     ]
#     resp = requests.put(url=f"{fus_url_api}/resource/{resource_type_name}/id/{project}/members")
#     print(resp)
#
# for project in open_projects:
#     resp = requests.post(url=f'{fus_url_api}/resource/{resource_type_name}/id/{project}',
#                          headers=auth_headers)
#     print(resp)
#
#     member_policy = [
#         {
#             "Add or update access": {
#                 "value": {
#                     "member_type": "default_user",
#                     "access_level": open_resource_policy_name
#                 }
#             }
#         }
#     ]
#     resp = requests.put(url=f"{fus_url_api}/resource/{resource_type_name}/id/{project}/members")
#     print(resp)




# fus_client = hca.auth.AuthClient(swagger_url=fus_url)
# # create our resource, with applicable actions
# fus_client.post_v1_resources(resource_type_name=resource_type_name,actions=['dss:GetBundle', 'dss:GetFile'])
# for projects in protected_projects:
#     fus_client.post_v1_resource_id(resource_type_name=resource_type_name,resource_id=projects)
#     fus_client.put_v1_resource_id_memeber(resource_type_name=resource_type_name,resource_id=projects) # unable to pass in the member into the request,
#
    #need to understand the relation of actions // policy for a given resource_type