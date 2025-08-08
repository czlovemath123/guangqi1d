module io_out
use datastructure
use communication
use mathlib
use phylib
use passivescalars
use eos
use HDF5
use mpi
implicit none

character(len=32) :: dsetname(100),filename_head
integer :: output_var_number,ncell_var,nxsurf_var,nysurf_car
character(len=32) :: attr_variables(5)

contains

subroutine simulation_config()
    character(len=24) :: fmt1,fmt2,fmt3
    attr_variables=[character(len=32) :: 'Time','nblocks','blk_size_nx','blk_size_ny','nblocks_tree']
    fmt1='(A20,1I15)'
    fmt2='(A20,1ES15.5E2)'
    fmt3='(A20,L15)'
    if (rank==0) then
        call chdir(path_out)
        open(unit=14,file=trim(path_root)//'/configuration.dat',status='replace',action='write')
        call print_config(14)
        close(14)
        call chdir(path_root)
    end if
end subroutine simulation_config

subroutine print_config(i)
    integer :: i
    character(len=24) :: fmt1,fmt2,fmt3
    fmt1='(A20,1I15)'
    fmt2='(A20,1ES15.5E2)'
    fmt3='(A20,L15)'
    write(unit=i,fmt=*) '&config_nml'
    write(unit=i,fmt=fmt2) 'CFL_number=', CFL
    write(unit=i,fmt=fmt1) 'nd=', nd
    write(unit=i,fmt=fmt1) 'number_of_frame=', nframe
    write(unit=i,fmt=fmt2) 'coord1min=', n_domain(1)
    write(unit=i,fmt=fmt2) 'coord1max=', n_domain(2)
    write(unit=i,fmt=fmt2) 'coord2min=', n_domain(3)
    write(unit=i,fmt=fmt2) 'coord2max=', n_domain(4)
    write(unit=i,fmt=fmt1) 'ncoord1=', nx
    write(unit=i,fmt=fmt1) 'ncoord2=', ny
    write(unit=i,fmt=fmt1) 'ncoord3=', nz
    write(unit=i,fmt=fmt1) 'igravity=', igravity
    write(unit=i,fmt=fmt1) 'igeometry=', igeometry
    write(unit=i,fmt=fmt1) 'icooling=', icooling
    write(unit=i,fmt=fmt1) 'iradiation=', iradiation
    write(unit=i,fmt=fmt1) 'maximum logical level', llevel_max
    write(unit=i,fmt=fmt2) 'mean atomic weight=', maw
    write(unit=i,fmt=fmt2) 'gamma_gas=', gamma_gas 
    write(unit=i,fmt=*) '/'
end subroutine print_config

subroutine output_blocks_parallel()
    character(len=32) :: filename_hdf,filename_xml
    character(len=16) :: fmt1,s1
    integer :: i,ierr,comm,info,error
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL
    fmt1='(I5.5)'
    write(s1,fmt1) iframe
    filename_hdf=trim(filename_head)//trim(s1)//'.h5'
    filename_xml=trim(filename_head)//trim(s1)//'.xdmf'
    call calculate_output_data()
    call hdf_output_block_parallel(filename_hdf)
    if (rank==0) then
        open(unit=89,file="last.dat",status="replace",action="write")
        write(89,fmt='(I6)') iframe
        close(89)
    end if
end subroutine output_blocks_parallel

subroutine calculate_output_data()
    !calculate some derivative data here
    type(blockdef), pointer :: blk
    integer :: i
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        if (iradiation/=0) then
            call calculate_entropy_block(blk)
        end if
        blk=>blk%next
    end do
end subroutine calculate_output_data

subroutine dsetname_filename()
    character(len=32) :: filename2
    integer :: i
    call chdir(path_out)
    open(unit=14,file=trim(path_root)//'/output_var_info.dat',status='old',action='read')
    read(unit=14,fmt=*) output_var_number
    do i=1,output_var_number
        read(unit=14,fmt=*) dsetname(i)
    end do
    if (iradiation==1) then
    else if (iradiation==4.or.iradiation==3) then
        dsetname(output_var_number+1)='Erad'
        dsetname(output_var_number+2)='Fradx'
        dsetname(output_var_number+3)='Rosseland'
        dsetname(output_var_number+4)='Planck'
        dsetname(output_var_number+5)='entropy'
        output_var_number=output_var_number+5
        if (larad) then
            dsetname(output_var_number+1)='aradx'
            output_var_number=output_var_number+1
        end if
    end if
#if     ieos==2
#if     ieosmodule==1
    dsetname(output_var_number+1)='H2'
    dsetname(output_var_number+2)='HI'
    dsetname(output_var_number+3)='HII'
    output_var_number=output_var_number+3
#elif   ieosmodule==2
    dsetname(output_var_number+1)='HI'
    dsetname(output_var_number+2)='HII'
    output_var_number=output_var_number+2
#elif   ieosmodule==3
    dsetname(output_var_number+1)='H2'
    dsetname(output_var_number+2)='HI'
    dsetname(output_var_number+3)='HII'
    dsetname(output_var_number+4)='HeI'
    dsetname(output_var_number+5)='HeII'
    dsetname(output_var_number+6)='HeIII'
    output_var_number=output_var_number+6
#elif   ieosmodule==4
    dsetname(output_var_number+1)='H2'
    dsetname(output_var_number+2)='HI'
    output_var_number=output_var_number+2
#endif
#endif
    if (lam_con) then
        if (igeometry==2) then
            dsetname(output_var_number+1)='omega'
            output_var_number=output_var_number+1
        else if(igeometry==3) then
            dsetname(output_var_number+1)='vphi'
            dsetname(output_var_number+2)='omega'
            output_var_number=output_var_number+2
        end if
    end if
    read(unit=14,fmt=*) filename2
    close(14)
    filename_head=trim(filename2)
    call chdir(path_root)
end subroutine dsetname_filename

subroutine hdf_output_block_parallel(filename)
    type(blockdef), pointer :: blk
    integer :: i,j,k,ijk(3),ip,comm,info
    integer :: error,dset_rank,arank
    integer(HSIZE_T), dimension(:), allocatable :: dims,dims_xsurface,dims_ysurface,dims_mg
    integer(HSIZE_T) :: blk_xdim(1),blk_ydim(1),att_dim(1)
    real(8) :: rho,temp,p,p_temp,log10_temp,log10_egv,eg,cs,xyz(3),rthetaphi(3),egv,nhtot,nhetot,species4(4),species5(5)
    character(len=16) :: fmt1,s1,passivename
    character(len=32) :: filename,groupname,blockname,attname
    character(len=128) :: alert,att_str
    integer(HID_T) :: file_id,group_id,dset_id,dspace_id,aspace_id,attr_id,plist_id1,plist_id2
    call chdir(path_out)
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL
    call h5open_f(error)
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id1, error)
    CALL h5pset_fapl_mpio_f(plist_id1, comm, info, error)
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error,access_prp=plist_id1)
    call h5pclose_f(plist_id1,error)
    dset_rank=2
    allocate(dims(dset_rank),dims_xsurface(dset_rank),dims_ysurface(dset_rank))
    dims=(/blk_size_nx,blk_size_ny/)
    blk_xdim=blk_size_nx+1
    dims_xsurface(1)=blk_size_nx+1
    blk_ydim=1
    dims_xsurface(2)=1
    arank=1
    att_dim(1)=1
    call blk_traversal(eos_species)
    call create_simulation_attr_parallel(file_id)
    call write_simulation_attr_parallel(file_id)
    do k=1,nblk_total
        fmt1='(I5.5)'
        write(s1,fmt1) k
        groupname='blk'//trim(s1)
        call h5gcreate_f(file_id,groupname,group_id,error)
        !!create the datasets
        call h5screate_simple_f(1,att_dim,dspace_id,error)
        call h5dcreate_f(group_id,'cpu_id',h5t_native_integer,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        call h5screate_simple_f(1,blk_xdim,dspace_id,error)
        call h5dcreate_f(group_id,'mesh_x',h5t_native_double,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        call h5screate_simple_f(1,blk_xdim-1,dspace_id,error)
        call h5dcreate_f(group_id,'x_center',h5t_native_double,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        do i=1,output_var_number
            call create_blk_dset_cell_parallel(blk,group_id,dims,dsetname(i))
            call create_blk_dset_xsurface_parallel(blk,group_id,dims_xsurface,dsetname(i))
            call create_blk_dset_ysurface_parallel(blk,group_id,dims_ysurface,dsetname(i))
        end do
        call h5gclose_f(group_id,error)
    end do
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id2,error)
    call h5pset_dxpl_mpio_f(plist_id2,h5fd_mpio_independent_f,error)
    blk=>llist_head
    do k=1,np_nblk(rank+1)
        fmt1='(I5.5)'
        write(s1,fmt1) blk%blk_id
        groupname='blk'//trim(s1)
        call h5gopen_f(file_id,groupname,group_id,error)
        call h5dopen_f(group_id,'cpu_id',dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_integer,rank,att_dim,error,xfer_prp=plist_id2)
        call h5dclose_f(dset_id,error)
        call h5dopen_f(group_id,'mesh_x',dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_double,blk%mesh_x,blk_xdim,error,xfer_prp=plist_id2)
        call h5dclose_f(dset_id,error)
        call h5dopen_f(group_id,'x_center',dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_double,blk%x_center(1:blk_size_nx),blk_xdim-1,error,xfer_prp=plist_id2)
        call h5dclose_f(dset_id,error)
        do i=1,output_var_number
            call write_blk_dset_cell_parallel(blk,group_id,dims,dsetname(i),plist_id2)
            call write_blk_dset_xsurface_parallel(blk,group_id,dims_xsurface,dsetname(i),plist_id2)
            call write_blk_dset_ysurface_parallel(blk,group_id,dims_ysurface,dsetname(i),plist_id2)
        end do
        call h5gclose_f(group_id,error)
        blk=>blk%next
    end do
    deallocate(dims,dims_xsurface)
    call h5fclose_f(file_id,error)
    call h5pclose_f(plist_id2, error)
    call h5close_f(error)
    call chdir(path_root)
end subroutine hdf_output_block_parallel

subroutine create_blk_dset_xsurface_parallel(blk,group_id,dims,dsname)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: xsurf_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id
    character(len=32) :: dsname
    integer :: error,ip,loc
    loc=findloc(xsurf_variables,dsname,dim=1)
    if (loc/=0) then
        call h5screate_simple_f(2,dims,dspace_id,error)
        call h5dcreate_f(group_id,trim(dsname),h5t_native_real,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        call h5sclose_f(dspace_id,error)
    end if
end subroutine create_blk_dset_xsurface_parallel

subroutine write_blk_dset_xsurface_parallel(blk,group_id,dims,dsname,plist_id)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: xsurf_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id,plist_id
    character(len=32) :: dsname
    integer :: error,ip
    logical :: xsurf
    allocate(xsurf_data(0:blk_size_nx,blk_size_ny))
    xsurf=.true.
    select case(dsname)
    case ('Fradx')
        xsurf_data=blk%Fradx(0:blk_size_nx,1:blk_size_ny,1)
    case ('aradx')
        xsurf_data=blk%aradx(0:blk_size_nx,1:blk_size_ny,1)
    case ('viscous_heat_xz')
        xsurf_data=blk%viscous_heat_xz(0:blk_size_nx,1:blk_size_ny,1)
    case ('torque_xz')
        xsurf_data=blk%torque_xz(0:blk_size_nx,1:blk_size_ny,1)
    case ('mass_xflux')
        xsurf_data=blk%xflux(1,0:blk_size_nx,1:blk_size_ny,1)
    case ('am_xflux')
        ip=ipassive(pn_am)
        xsurf_data=blk%xpflux(ip,0:blk_size_nx,1:blk_size_ny,1)
    case ('torque_work')
    case ('hydro_energy_xflux')
        xsurf_data=blk%xflux(5,0:blk_size_nx,1:blk_size_ny,1)
    case default
        xsurf=.false.
    end select
    if (xsurf) then
        call h5dopen_f(group_id,trim(dsname),dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_real,xsurf_data,dims,error,xfer_prp=plist_id)
        call h5dclose_f(dset_id,error)
    end if
    deallocate(xsurf_data)
end subroutine write_blk_dset_xsurface_parallel

subroutine create_blk_dset_ysurface_parallel(blk,group_id,dims,dsname)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: ysurf_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id
    character(len=32) :: dsname
    integer :: error,loc
    loc=findloc(ysurf_variables,dsname,dim=1)
    if (loc/=0) then
        call h5screate_simple_f(2,dims,dspace_id,error)
        call h5dcreate_f(group_id,trim(dsname),h5t_native_real,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        call h5sclose_f(dspace_id,error)
    end if
end subroutine create_blk_dset_ysurface_parallel

subroutine write_blk_dset_ysurface_parallel(blk,group_id,dims,dsname,plist_id)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: ysurf_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id,plist_id
    character(len=32) :: dsname
    integer :: error
    logical :: ysurf
    allocate(ysurf_data(blk_size_nx,0:blk_size_ny))
    ysurf=.true.
    select case(dsname)
    case ('Frady')
        ysurf_data=blk%Frady(1:blk_size_nx,0:blk_size_ny,1)
    case ('viscous_heat_yz')
        ysurf_data=blk%viscous_heat_yz(1:blk_size_nx,0:blk_size_ny,1)
    case ('torque_yz')
        ysurf_data=blk%torque_yz(1:blk_size_nx,0:blk_size_ny,1)
    case ('mass_yflux')
        ysurf_data=blk%yflux(1,1:blk_size_nx,0:blk_size_ny,1)
    case ('am_yflux')
        ysurf_data=blk%ypflux(1,1:blk_size_nx,0:blk_size_ny,1)
    case default
        ysurf=.false.
    end select
    if (ysurf) then
        call h5dopen_f(group_id,trim(dsname),dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_real,ysurf_data,dims,error,xfer_prp=plist_id)
        call h5dclose_f(dset_id,error)
    end if
    deallocate(ysurf_data)
end subroutine write_blk_dset_ysurface_parallel

subroutine create_blk_dset_cell_parallel(blk,group_id,dims,dsname)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: cell_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id
    character(len=32) :: dsname
    character(len=16) :: passivename
    integer :: ip,error,loc
    loc=findloc(cell_variables,dsname,dim=1)
    if (loc/=0) then
        call h5screate_simple_f(2,dims,dspace_id,error)
        call h5dcreate_f(group_id,trim(dsname),h5t_native_real,dspace_id,dset_id,error)
        call h5dclose_f(dset_id,error)
        call h5sclose_f(dspace_id,error)
    end if
end subroutine create_blk_dset_cell_parallel

subroutine write_blk_dset_cell_parallel(blk,group_id,dims,dsname,plist_id)
    type(blockdef), pointer :: blk
    real(4), dimension(:,:), allocatable :: cell_data
    integer(HSIZE_T), dimension(:), allocatable :: dims
    integer(HID_T) :: group_id,dspace_id,dset_id,plist_id
    character(len=32) :: dsname
    character(len=16) :: passivename
    integer :: ip,error
    logical :: cell
    allocate(cell_data(blk_size_nx,blk_size_ny))
    cell=.true.
    select case(dsname)
    case ('rho')
        cell_data=blk%w(1,1:blk_size_nx,1:blk_size_ny,1)
    case ('vx')
        cell_data=blk%w(2,1:blk_size_nx,1:blk_size_ny,1)
    case ('vy')
        cell_data=blk%w(3,1:blk_size_nx,1:blk_size_ny,1)
    case ('vz')
        cell_data=blk%w(4,1:blk_size_nx,1:blk_size_ny,1)
    case ('pres')
        cell_data=blk%w(5,1:blk_size_nx,1:blk_size_ny,1)
    case ('cs')
        cell_data=blk%cs(1:blk_size_nx,1:blk_size_ny,1)
    case ('temp')
        cell_data=blk%temp(1:blk_size_nx,1:blk_size_ny,1)
    case ('egv')
        cell_data=blk%egv(1:blk_size_nx,1:blk_size_ny,1)
    case ('entropy')
        cell_data=blk%entropy(1:blk_size_nx,1:blk_size_ny,1)
    case ('vphi')
        passivename='vphi'
        ip=ipassive(passivename)
        cell_data=blk%passive_scalar(ip,1:blk_size_nx,1:blk_size_ny,1)
    case ('omega')
        cell_data=blk%omega(1:blk_size_nx,1:blk_size_ny,1)
    case ('irrad')
        cell_data=blk%irrad(1:blk_size_nx,1:blk_size_ny,1)
    case ('H2')
        cell_data=blk%H2(1:blk_size_nx,1:blk_size_ny,1)
    case ('HI')
        cell_data=blk%HI(1:blk_size_nx,1:blk_size_ny,1)
    case ('HII')
        cell_data=blk%HII(1:blk_size_nx,1:blk_size_ny,1)
    case ('HeI')
        cell_data=blk%HeI(1:blk_size_nx,1:blk_size_ny,1)
    case ('HeII')
        cell_data=blk%HeII(1:blk_size_nx,1:blk_size_ny,1)
    case ('HeIII')
        cell_data=blk%HeIII(1:blk_size_nx,1:blk_size_ny,1)
    case ('Rosseland')
        cell_data=blk%kappa_rosseland(1:blk_size_nx,1:blk_size_ny,1)
    case ('sigma_Rosseland')
        cell_data=blk%sigma_rosseland(1:blk_size_nx,1:blk_size_ny,1)
    case ('Planck')
        cell_data=blk%kappa_planck(1:blk_size_nx,1:blk_size_ny,1)
    case ('sigma_Planck')
        cell_data=blk%sigma_planck(1:blk_size_nx,1:blk_size_ny,1)
    case ('viscous_heat')
        cell_data=blk%viscous_heat(1:blk_size_nx,1:blk_size_ny,1)
    case ('Erad')
        cell_data=blk%Erad(1:blk_size_nx,1:blk_size_ny,1)
    case ('Erad_int')
        cell_data=blk%Erad_int(1:blk_size_nx,1:blk_size_ny,1)
    case ('Erad_source')
        cell_data=blk%Erad_source(1:blk_size_nx,1:blk_size_ny,1)
    case ('isotemp')
        cell_data=blk%iso_temp(1:blk_size_nx,1:blk_size_ny,1)
    case ('fld_heating')
        cell_data=blk%fld_heating(1:blk_size_nx,1:blk_size_ny,1)
    case ('fld_heating_viscous')
        cell_data=blk%fld_heating_viscous(1:blk_size_nx,1:blk_size_ny,1)
    case default
        cell=.false.
    end select
    if (cell) then
        call h5dopen_f(group_id,trim(dsname),dset_id,error)
        call h5dwrite_f(dset_id,h5t_native_real,cell_data,dims,error,xfer_prp=plist_id)
        call h5dclose_f(dset_id,error)
    end if
    deallocate(cell_data)
end subroutine write_blk_dset_cell_parallel

subroutine create_simulation_attr_double_parallel(file_id,str,val)
    integer(HID_T) :: file_id,aspace_id,attr_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: arank,error
    real(8) :: val
    character(len=*) :: str
    arank=1
    att_dim(1)=1
    call h5screate_simple_f(arank,att_dim,aspace_id,error)
    call h5acreate_f(file_id,trim(str),h5t_native_double,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
end subroutine create_simulation_attr_double_parallel

subroutine create_simulation_attr_int_parallel(file_id,str,val)
    integer(HID_T) :: file_id,aspace_id,attr_id
    integer(HSIZE_T) :: att_dim(1)
    integer :: arank,error,val
    character(len=*) :: str
    arank=1
    att_dim(1)=1
    call h5screate_simple_f(arank,att_dim,aspace_id,error)
    call h5acreate_f(file_id,trim(str),h5t_native_integer,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
end subroutine create_simulation_attr_int_parallel

subroutine create_simulation_attr_parallel(file_id)
    integer(HID_T) :: file_id,aspace_id,attr_id,dspace_id,dset_id
    integer(HSIZE_T) :: att_dim(1),dims(2),dim_level(1)
    integer :: i,arank,error,val
    att_dim(1)=1
    call h5screate_simple_f(1,att_dim,aspace_id,error)
    call h5acreate_f(file_id,'Time',h5t_native_double,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
    call h5screate_simple_f(1,att_dim,aspace_id,error)
    call h5acreate_f(file_id,'nblocks',h5t_native_integer,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
    call h5screate_simple_f(1,att_dim,aspace_id,error)
    call h5acreate_f(file_id,'blk_size_nx',h5t_native_integer,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
    call h5screate_simple_f(1,att_dim,aspace_id,error)
    call h5acreate_f(file_id,'blk_size_ny',h5t_native_integer,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
    call h5screate_simple_f(1,att_dim,aspace_id,error)
    call h5acreate_f(file_id,'nblocks_tree',h5t_native_integer,aspace_id,attr_id,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id, error)
    nblocks_tree=1
    call count_child(blk_root,nblocks_tree)
    dims=(/nblocks_tree,6/)
    call h5screate_simple_f(2,dims,dspace_id,error)
    call h5dcreate_f(file_id,'block_tree',h5t_native_integer,dspace_id,dset_id,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)
    dim_level(1)=nblk_total
    call h5screate_simple_f(1,dim_level,dspace_id,error)
    call h5dcreate_f(file_id,'levels',h5t_native_integer,dspace_id,dset_id,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)
end subroutine create_simulation_attr_parallel

subroutine write_simulation_attr_parallel(file_id)
    type(blockdef), pointer :: blk
    integer(HID_T) :: file_id,aspace_id,attr_id,dset_id
    integer(HSIZE_T) :: att_dim(1),dims(2),dim_level(1)
    integer :: i,arank,error,val,istart,iend,req,stat(MPI_STATUS_SIZE),ierr
    integer, dimension(:,:), allocatable :: keys
    integer, dimension(:), allocatable :: level_local,level_all
    att_dim(1)=1
    nblocks_tree=1
    call count_child(blk_root,nblocks_tree)
    allocate(keys(nblocks_tree,6),level_local(np_nblk(rank+1)),level_all(nblk_total))
    keys=0
    dims=(/nblocks_tree,6/)
    call save_the_keys(keys)
    call h5dopen_f(file_id,'block_tree',dset_id,error)
    call h5dwrite_f(dset_id,h5t_native_integer,keys,dims,error)
    call h5dclose_f(dset_id,error)
    call h5aopen_f(file_id,'Time',attr_id,error)
    call h5awrite_f(attr_id,h5t_native_double,time_sys%t,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5aopen_f(file_id,'nblocks',attr_id,error)
    call h5awrite_f(attr_id,h5t_native_integer,nblk_total,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5aopen_f(file_id,'blk_size_nx',attr_id,error)
    call h5awrite_f(attr_id,h5t_native_integer,blk_size_nx,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5aopen_f(file_id,'blk_size_ny',attr_id,error)
    call h5awrite_f(attr_id,h5t_native_integer,blk_size_ny,att_dim,error)
    call h5aclose_f(attr_id,error)
    call h5aopen_f(file_id,'nblocks_tree',attr_id,error)
    call h5awrite_f(attr_id,h5t_native_integer,nblocks_tree,att_dim,error)
    call h5aclose_f(attr_id,error)
    blk=>llist_head
    do i=1,np_nblk(rank+1)
        level_local(i)=blk%level
        blk=>blk%next
    end do
    if (np>1) then
        if (rank/=0) then
            call mpi_isend(level_local,np_nblk(rank+1),mpi_integer,0,1,MPI_COMM_WORLD,req,ierr)
        else
            level_all(1:np_nblk(1))=level_local
            do i=2,np
                istart=sum(np_nblk(1:i-1))+1
                iend=istart+np_nblk(i)-1
                call mpi_recv(level_all(istart:iend),np_nblk(i),mpi_integer,i-1,1,MPI_COMM_WORLD,stat,ierr)
            end do
        end if
    else
        level_all=level_local
    end if
    call mpi_bcast(level_all,nblk_total,mpi_integer,0,MPI_COMM_WORLD,ierr)
    call h5dopen_f(file_id,'levels',dset_id,error)
    dim_level(1)=nblk_total
    call h5dwrite_f(dset_id,h5t_native_integer,level_all,dim_level,error)
    call h5dclose_f(dset_id,error)
    deallocate(keys,level_local,level_all)
end subroutine write_simulation_attr_parallel

recursive subroutine count_child(blk,n)
    !preorder traversal
    type(blockdef), pointer :: blk
    integer :: idx,idy,n
    idy=1
    do idx=1,2
        if (associated(blk%children(idx,idy)%blk)) then
            n=n+1;call count_child(blk%children(idx,idy)%blk,n)
        end if
    end do
end subroutine count_child

subroutine save_the_keys(keys)
    integer, dimension(:,:), allocatable :: keys
    integer :: i
    i=1
    call save_key(blk_root,keys,i)
end subroutine save_the_keys

recursive subroutine save_key(blk,keys,i)
    !preorder traversal
    type(blockdef), pointer :: blk
    integer, dimension(:,:), allocatable :: keys
    integer :: i,key_p(3),key(3),idx,idy
    if (associated(blk%pa)) then
        key_p=blk%pa%key
    else
        key_p=blk%key
    end if
    key=blk%key
    keys(i,1:6)=(/key_p(1:3),key(1:3)/)
    idy=1
    do idx=1,2
        if (associated(blk%children(idx,idy)%blk)) then
            i=i+1;call save_key(blk%children(idx,idy)%blk,keys,i)
        end if
    end do
end subroutine save_key

subroutine create_xdmf_xml(filename_hdf,filename_xml)
    character(len=32) :: fmt1,s1,filename_hdf,filename_xml,blockname
    type(blockdef), pointer :: blk
    integer :: i,j,k
    call chdir(path_out)
    fmt1='(I5.5)'
    open(unit=15,file=filename_xml,status='replace',action='write')
    write(15,'(a)')'<Xdmf Version="2.0">'
    write(15,'(a)')'<Information Name="TimeVaryingMetaData" Value="True"/>'
    write(15,'(a)')'<Domain Name="domain">'
    write(15,'(a)')'<Grid Name="mesh" GridType="Collection">'
    write(15,*)' <Time Value="',time_sys%t,'"/>'
    blk=>blk_head
    do j=1,nblk_total
        write(s1,fmt1) blk%blk_id
        blockname='blk'//trim(s1)
        write(15,*)'    <Grid Name="',trim(blockname),'" GridType="Uniform">'
        write(15,*)'    <Topology TopologyType="2DRectMesh" NumberOfElements="',blk_size_ny+1,blk_size_nx+1,'"/>'
        write(15,*)'    <Geometry GeometryType="VXVY">'
        write(15,*)'        <DataItem Dimensions="',blk_size_nx+1,'" NumberType="Double" Precision="8" Format="HDF">'
        write(15,*)'            '//trim(filename_hdf)//':/'//trim(blockname)//'/'//'mesh_x'
        write(15,*)'        </DataItem>'
        write(15,*)'        <DataItem Dimensions="',blk_size_ny+1,'" NumberType="Double" Precision="8" Format="HDF">'
        write(15,*)'            '//trim(filename_hdf)//':/'//trim(blockname)//'/'//'mesh_y'
        write(15,*)'        </DataItem>'
        write(15,*)'    </Geometry>'
        do i=1,output_var_number
            write(15,*)'    <Attribute Name="',trim(dsetname(i)),'" AttributeType="Scalar" Center="Cell">'
            write(15,*)'        <DataItem Dimensions="',blk_size_ny,blk_size_nx,'" NumberType="Double" Precision="8" Format="HDF">'
            write(15,*)'        '//trim(filename_hdf)//':/'//trim(blockname)//'/'//trim(dsetname(i))
            write(15,*)'        </DataItem>'
            write(15,*)'    </Attribute>'
        end do
        write(15,*)'    </Grid>'
        blk=>blk%next
    end do
    write(15,'(a)')'</Grid>'
    write(15,'(a)')'</Domain>'
    write(15,'(a)')'</Xdmf>'
    close(15)
    call chdir(path_root)
end subroutine create_xdmf_xml

subroutine status_quo_columns()
    if (refine_type=='adaptive'.or.refine_type=='mixed') then
        if (rank==0) write(*,'(A8,5A18)')'frame','next output time','t_final','frame dt','nblks/np','load ratio'
    else
        if (rank==0) write(*,'(A8,4A18)')'frame','next output time','t_final','frame dt','nblks/np'
    end if
end subroutine status_quo_columns

subroutine status_quo(t1,t2)
    real(8) :: t1,t2
    if (refine_type=='adaptive'.or.refine_type=='mixed') then
        if (rank==0) then
            write(*,'(I8,3ES18.6E2,2F18.2)') iframe,time_sys%t_next/timescale,time_sys%t_final/timescale,t2-t1,real(nblk_total)/real(np), &
            real(maxval(np_nblk))/real(minval(np_nblk))
        end if
    else
        if (rank==0) then
            write(*,'(I8,3ES18.6E2,F18.2)') iframe,time_sys%t_next/timescale,time_sys%t_final/timescale,t2-t1,real(nblk_total)/real(np)
        end if
    end if
end subroutine status_quo

end module io_out
