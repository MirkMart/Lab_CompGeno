# Server_Didattica

This folder aims to help managing and organise Didattica's server.

SSh connection:

```bash
# X is the number if present
ssh PERSONALE^nome.cognomeX@137.204.142.242
```

Server's characteristics:

- SO: Ubuntu desktop 24.04.2 LTS (Noble Numbat).
- Architecture: x86_64 (64-bit intel/AMD)
- Storage:
  - Filesystem      Size   Mounted_on
  - tmpfs            76G   /run
  - /dev/sdh1       440G   /
  - tmpfs           378G   /dev/shm
  - tmpfs           5.0M   /run/lock
  - /dev/sdg1       916G   /Store1Tb
  - /dev/md0        2.7T   /home
  - tmpfs           92K    /run/user/120
  - tmpfs           80K    /run/user/10005155
  - tmpfs           96K    /run/user/10997013
