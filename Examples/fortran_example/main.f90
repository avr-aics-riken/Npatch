
!/*
! * Nagata Patch Library
! *
! *
! * Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
! * All rights reserved.
! *
! *
! */

!////////////////////////////////////////////////////////////////////////////
!///
!/// Fortranインターフェース サンプル
!///   三角形の頂点列データを入力として、以下の処理を行う
!///     ・頂点の法線ベクトルを求める
!///     ・長田パッチの生成
!///     ・三角形の各辺の中点の長田パッチ上の曲面補間点を求める
!///     ・三角形の頂点と上記の曲面補間点より、三角形数を４倍としたデータを生成する
!///   - 当サンプルは頂点の法線ベクトルを求めるのに単純にループさせているが
!///     三角形数が多い場合には処理時間がかかるようになるので注意
!///
!////////////////////////////////////////////////////////////////////////////


! 精度の定義
!   gfortranの場合、includeファイル内ではプリプロセッサが効かず
!   コンパイルエラーとなるため、それぞれのファイルで#defineする
#ifdef _REAL_IS_DOUBLE_
#define NPT_REAL_PN   8
#else
#define NPT_REAL_PN   4
#endif

#define NMAX 20
#define NPT_ALW_L 0.001

! 外積
subroutine FCalcOutProduct( vec1, vec2, vec_o )
    implicit none
    real(NPT_REAL_PN), intent(in)  :: vec1(3)
    real(NPT_REAL_PN), intent(in)  :: vec2(3)
    real(NPT_REAL_PN), intent(out) :: vec_o(3)

    vec_o(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    vec_o(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    vec_o(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    return
end subroutine FCalcOutProduct

! ベクトルの正規化
subroutine FCalcNormalize( vec )
    implicit none
    real(NPT_REAL_PN), intent(inout)  :: vec(3)
    real(NPT_REAL_PN) :: len
    len = sqrt( vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3) )
    vec(1) = vec(1)/len
    vec(2) = vec(2)/len
    vec(3) = vec(3)/len
end subroutine FCalcNormalize

! 面の法線ベクトル
subroutine get_tri_normal( p1, p2, p3, norm )
    implicit none
    real(NPT_REAL_PN), intent(in)  :: p1(3)
    real(NPT_REAL_PN), intent(in)  :: p2(3)
    real(NPT_REAL_PN), intent(in)  :: p3(3)
    real(NPT_REAL_PN), intent(out) :: norm(3)
    real(NPT_REAL_PN) :: vec12(3), vec13(3);

    vec12(1)=p2(1)-p1(1);  vec12(2)=p2(2)-p1(2); vec12(3)=p2(3)-p1(3);
    vec13(1)=p3(1)-p1(1);  vec13(2)=p3(2)-p1(2); vec13(3)=p3(3)-p1(3);
    call FCalcOutProduct( vec12, vec13, norm );
    call FCalcNormalize( norm );
end subroutine get_tri_normal

! STLファイル出力
subroutine output_stl_file( file_name, num, tri, norm )
    implicit none
    character(64), intent(in)  :: file_name
    integer, intent(in)  :: num
    real(NPT_REAL_PN), intent(in) :: tri(3,3,num)
    real(NPT_REAL_PN), intent(in) :: norm(3,num)

    integer :: iunit = 7
    integer :: i,j

    open( iunit, file=file_name )

    write(iunit,'("solid model1")')
    do i=1,num
        write(iunit,'("    facet normal ",f10.6,1x,f10.6,1x,f10.6)') norm(1,i),norm(2,i),norm(3,i)
        write(iunit,'("    outer loop")');
        do j=1,3
            write(iunit,'("        vertex ",f12.5,1x,f12.5,1x,f12.5)') tri(1,j,i),tri(2,j,i),tri(3,j,i)
        enddo
        write(iunit,'("    endloop")')
        write(iunit,'("    endfacet")')
    enddo
    write(iunit,'("endsolid model1")')

    close(iunit)
end subroutine output_stl_file

! NPTファイル出力
subroutine output_npt_file( file_name, num, tri, npatch )
    implicit none
    character(64), intent(in)  :: file_name
    integer, intent(in)  :: num
    real(NPT_REAL_PN), intent(in) :: tri   (3,3,num)
    real(NPT_REAL_PN), intent(in) :: npatch(3,7,num)

    integer :: iunit = 8
    integer :: i,j

    open( iunit, file=file_name )

    write(iunit,'(i10)') num
    do i=1,num
        write(iunit,'("facet")')
        do j=1,3
            write(iunit,'("    vertex ",f12.5,1x,f12.5,1x,f12.5)') tri(1,j,i),tri(2,j,i),tri(3,j,i)
        enddo
        do j=1,7
            write(iunit,'("    coef",i1,f12.5,1x,f12.5,1x,f12.5)') j,npatch(1,j,i),npatch(2,j,i),npatch(3,j,i)
        enddo
    enddo

    close(iunit)
end subroutine output_npt_file


!//----------------------------------------------------
!//  メインルーチン
!//----------------------------------------------------

program main
    implicit none
    real(NPT_REAL_PN) :: p0(3),p1(3),p2(3),p3(3),p4(3),p5(3),p6(3),p7(3)
#if 1
    real(NPT_REAL_PN) :: midpnt(3)  ! 中間点
#endif
    integer :: num_tri;
    real(NPT_REAL_PN) :: tri       (3,3,NMAX)  ! 三角形の座標
    real(NPT_REAL_PN) :: plane_norm(3,NMAX)    ! 三角形の法線ベクトル
    real(NPT_REAL_PN) :: vtx_norm  (3,3,NMAX)  ! 三角形の頂点の法線ベクトル
    real(NPT_REAL_PN) :: npatch    (3,7,NMAX)  ! 長田パッチパラメータ

    integer :: i,j, iret;
    integer :: ip,jp;
    real(NPT_REAL_PN) :: norm_tmp(3)
    real(NPT_REAL_PN) :: eps = 0.01;
    integer :: num_4
    real(NPT_REAL_PN) :: eta, xi;
    real(NPT_REAL_PN) :: p12(3), p23(3),p31(3)
    real(NPT_REAL_PN) :: tri_4      (3,3,4*NMAX)  ! 長田パッチ補間により４倍したポリゴン
    real(NPT_REAL_PN) :: plane_norm4(3,4*NMAX)    ! ４倍したポリゴンの法線ベクトル

    character(64) :: file_name_stl_in   ='stl_in.stl'
    character(64) :: file_name_stl4_out ='stl4_out.stl'
    character(64) :: file_name_npt_out  ='npt_out.npt'

    write(*,*) '#### Program Start'

    ! 頂点の設定
#if 0
    ! 三角形３個（開放端あり）
    num_tri = 3
    p0(1)=  0.0;  p0(2)=  0.0;  p0(3)=  0.0
    p1(1)=100.0;  p1(2)=  0.0;  p1(3)=  0.0
    p2(1)=  0.0;  p2(2)=100.0;  p2(3)=  0.0
    p3(1)=  0.0;  p3(2)=  0.0;  p3(3)=100.0

    tri(1,1,1)=p0(1);  tri(2,1,1)=p0(2);  tri(3,1,1)=p0(3)
    tri(1,2,1)=p2(1);  tri(2,2,1)=p2(2);  tri(3,2,1)=p2(3)
    tri(1,3,1)=p1(1);  tri(2,3,1)=p1(2);  tri(3,3,1)=p1(3)

    tri(1,1,2)=p0(1);  tri(2,1,2)=p0(2);  tri(3,1,2)=p0(3);
    tri(1,2,2)=p1(1);  tri(2,2,2)=p1(2);  tri(3,2,2)=p1(3);
    tri(1,3,2)=p3(1);  tri(2,3,2)=p3(2);  tri(3,3,2)=p3(3);

    tri(1,1,3)=p0(1);  tri(2,1,3)=p0(2);  tri(3,1,3)=p0(3);
    tri(1,2,3)=p3(1);  tri(2,2,3)=p3(2);  tri(3,2,3)=p3(3);
    tri(1,3,3)=p2(1);  tri(2,3,3)=p2(2);  tri(3,3,3)=p2(3);
#endif

#if 1
    ! 直方体
    num_tri = 12;
    p0(1)=  0.0;  p0(2)=  0.0;  p0(3)=  0.0
    p1(1)=100.0;  p1(2)=  0.0;  p1(3)=  0.0
    p2(1)=100.0;  p2(2)=100.0;  p2(3)=  0.0
    p3(1)=  0.0;  p3(2)=100.0;  p3(3)=  0.0
    p4(1)=  0.0;  p4(2)=  0.0;  p4(3)=100.0
    p5(1)=100.0;  p5(2)=  0.0;  p5(3)=100.0
    p6(1)=100.0;  p6(2)=100.0;  p6(3)=100.0
    p7(1)=  0.0;  p7(2)=100.0;  p7(3)=100.0
    ! 三角形の設定
    !     底面 p0->p3->p2->p1
    tri(1,1,1)=p0(1);  tri(2,1,1)=p0(2);  tri(3,1,1)=p0(3)
    tri(1,2,1)=p3(1);  tri(2,2,1)=p3(2);  tri(3,2,1)=p3(3)
    tri(1,3,1)=p1(1);  tri(2,3,1)=p1(2);  tri(3,3,1)=p1(3)

    tri(1,1,2)=p2(1);  tri(2,1,2)=p2(2);  tri(3,1,2)=p2(3)
    tri(1,2,2)=p1(1);  tri(2,2,2)=p1(2);  tri(3,2,2)=p1(3)
    tri(1,3,2)=p3(1);  tri(2,3,2)=p3(2);  tri(3,3,2)=p3(3)
    !     上面 p4->p5->p6->p7
    tri(1,1,3)=p4(1);  tri(2,1,3)=p4(2);  tri(3,1,3)=p4(3)
    tri(1,2,3)=p5(1);  tri(2,2,3)=p5(2);  tri(3,2,3)=p5(3)
    tri(1,3,3)=p7(1);  tri(2,3,3)=p7(2);  tri(3,3,3)=p7(3)

    tri(1,1,4)=p6(1);  tri(2,1,4)=p6(2);  tri(3,1,4)=p6(3)
    tri(1,2,4)=p7(1);  tri(2,2,4)=p7(2);  tri(3,2,4)=p7(3)
    tri(1,3,4)=p5(1);  tri(2,3,4)=p5(2);  tri(3,3,4)=p5(3)
    !     Front  p0->p1->p5->p4
    tri(1,1,5)=p0(1);  tri(2,1,5)=p0(2);  tri(3,1,5)=p0(3)
    tri(1,2,5)=p1(1);  tri(2,2,5)=p1(2);  tri(3,2,5)=p1(3)
    tri(1,3,5)=p4(1);  tri(2,3,5)=p4(2);  tri(3,3,5)=p4(3)

    tri(1,1,6)=p5(1);  tri(2,1,6)=p5(2);  tri(3,1,6)=p5(3)
    tri(1,2,6)=p4(1);  tri(2,2,6)=p4(2);  tri(3,2,6)=p4(3)
    tri(1,3,6)=p1(1);  tri(2,3,6)=p1(2);  tri(3,3,6)=p1(3)
    !     Right  p1->p2->p6->p5
    tri(1,1,7)=p1(1);  tri(2,1,7)=p1(2);  tri(3,1,7)=p1(3)
    tri(1,2,7)=p2(1);  tri(2,2,7)=p2(2);  tri(3,2,7)=p2(3)
    tri(1,3,7)=p5(1);  tri(2,3,7)=p5(2);  tri(3,3,7)=p5(3)

    tri(1,1,8)=p6(1);  tri(2,1,8)=p6(2);  tri(3,1,8)=p6(3)
    tri(1,2,8)=p5(1);  tri(2,2,8)=p5(2);  tri(3,2,8)=p5(3)
    tri(1,3,8)=p2(1);  tri(2,3,8)=p2(2);  tri(3,3,8)=p2(3)
    !     Back  p2->p3->p7->p6
    tri(1,1,9)=p2(1);  tri(2,1,9)=p2(2);  tri(3,1,9)=p2(3)
    tri(1,2,9)=p3(1);  tri(2,2,9)=p3(2);  tri(3,2,9)=p3(3)
    tri(1,3,9)=p6(1);  tri(2,3,9)=p6(2);  tri(3,3,9)=p6(3)

    tri(1,1,10)=p7(1);  tri(2,1,10)=p7(2);  tri(3,1,10)=p7(3)
    tri(1,2,10)=p6(1);  tri(2,2,10)=p6(2);  tri(3,2,10)=p6(3)
    tri(1,3,10)=p3(1);  tri(2,3,10)=p3(2);  tri(3,3,10)=p3(3)
    !     Left  p3->p0->p4->p7
    tri(1,1,11)=p3(1);  tri(2,1,11)=p3(2);  tri(3,1,11)=p3(3)
    tri(1,2,11)=p0(1);  tri(2,2,11)=p0(2);  tri(3,2,11)=p0(3)
    tri(1,3,11)=p7(1);  tri(2,3,11)=p7(2);  tri(3,3,11)=p7(3)

    tri(1,1,12)=p4(1);  tri(2,1,12)=p4(2);  tri(3,1,12)=p4(3)
    tri(1,2,12)=p7(1);  tri(2,2,12)=p7(2);  tri(3,2,12)=p7(3)
    tri(1,3,12)=p0(1);  tri(2,3,12)=p0(2);  tri(3,3,12)=p0(3)
#endif


    ! 面の法線ベクトル設定
    do i=1, num_tri
        call get_tri_normal( tri(1,1,i), tri(1,2,i), tri(1,3,i), plane_norm(1,i) )
    enddo

    ! 確認のためSTLファイルに出力
    call output_stl_file( file_name_stl_in, num_tri, tri, plane_norm )


    ! 頂点ベクトルの設定
    !    点数が少ないので単純ループ
    do i=1, num_tri  ! 3角形のloop
        do ip=1, 3      ! 頂点のloop
           norm_tmp(1) = plane_norm(1,i)   ! 自身の法線ベクトル設定
           norm_tmp(2) = plane_norm(2,i)
           norm_tmp(3) = plane_norm(3,i)

           do j=1, num_tri  ! 相手の3角形のloop
               if( i == j ) cycle;
               do jp=1, 3     ! 相手の頂点のloop

                   if( abs(tri(1,ip,i)-tri(1,jp,j)) < eps  .and. &
                       abs(tri(2,ip,i)-tri(2,jp,j)) < eps  .and. &
                       abs(tri(3,ip,i)-tri(3,jp,j)) < eps       ) then

                       norm_tmp(1) = norm_tmp(1) + plane_norm(1,j)   ! 接している面の法線ベクトル
                       norm_tmp(2) = norm_tmp(2) + plane_norm(2,j)
                       norm_tmp(3) = norm_tmp(3) + plane_norm(3,j)
                       exit
                   endif
               enddo
           enddo
           call FCalcNormalize( norm_tmp );
           vtx_norm(1,ip,i)=norm_tmp(1); vtx_norm(2,ip,i)=norm_tmp(2); vtx_norm(3,ip,i)=norm_tmp(3);
        enddo
    enddo

    ! 長田パッチ変換
    do i=1, num_tri    ! 3角形のloop
       call fnpt_param_crt (  &
                  tri(1,1,i), vtx_norm(1,1,i), &
                  tri(1,2,i), vtx_norm(1,2,i), &
                  tri(1,3,i), vtx_norm(1,3,i), &
                  npatch(1,1,i), npatch(1,2,i), npatch(1,3,i), npatch(1,4,i), &
                  npatch(1,5,i), npatch(1,6,i), npatch(1,7,i),  &
                  iret &
                )
       if( iret .ne. 1 ) then
           write(*,'("#### Error npt_param_crt() ret=",i3," i=",i9)') iret,i
           stop 1;
       endif
    enddo

    ! 確認のため長田パッチファイルに出力
    call output_npt_file( file_name_npt_out, num_tri, tri, npatch )

#if 0
    ! eta,xi取得テスト
    !    三角形１の辺１の中点
    midpnt(1)=( tri(1,1,1)+tri(1,2,1) )/2.0
    midpnt(2)=( tri(2,1,1)+tri(2,2,1) )/2.0
    midpnt(3)=( tri(3,1,1)+tri(3,2,1) )/2.0
    call fnpt_cvt_pos_to_eta_xi( midpnt, tri(1,1,1), tri(1,2,1), tri(1,3,1), eta, xi )
    write(*,'("---- edge1 mid point eta=",f10.6," xi=",f10.6)') eta,xi
    !    三角形１の辺２の中点
    midpnt(1)=( tri(1,2,1)+tri(1,3,1) )/2.0
    midpnt(2)=( tri(2,2,1)+tri(2,3,1) )/2.0
    midpnt(3)=( tri(3,2,1)+tri(3,3,1) )/2.0
    call fnpt_cvt_pos_to_eta_xi( midpnt, tri(1,1,1), tri(1,2,1), tri(1,3,1), eta, xi )
    write(*,'("---- edge2 mid point eta=",f10.6," xi=",f10.6)') eta,xi
    !    三角形１の辺３の中点
    midpnt(1)=( tri(1,3,1)+tri(1,1,1) )/2.0
    midpnt(2)=( tri(2,3,1)+tri(2,1,1) )/2.0
    midpnt(3)=( tri(3,3,1)+tri(3,1,1) )/2.0
    call fnpt_cvt_pos_to_eta_xi( midpnt, tri(1,1,1), tri(1,2,1), tri(1,3,1), eta, xi )
    write(*,'("---- edge3 mid point eta=",f10.6," xi=",f10.6)') eta,xi
#endif

    ! 各辺の中点を補間して点を接点として追加
    !    ポリゴン数が４倍になる
    num_4 = 0
    do i=1, num_tri    ! 3角形のloop
        ! 辺１の中点の曲面補間点
        eta=0.5; xi=0.0
        call fnpt_correct_pnt( &
                eta, xi,  &
                tri(1,1,i), tri(1,2,i), tri(1,3,i),  &
                npatch(1,1,i), npatch(1,2,i), npatch(1,3,i), npatch(1,4,i), &
                npatch(1,5,i), npatch(1,6,i), npatch(1,7,i),  &
                p12 &
             )

        ! 辺２の中点の曲面補間点
        eta=1.0; xi=0.5
        call fnpt_correct_pnt( &
                eta, xi,  &
                tri(1,1,i), tri(1,2,i), tri(1,3,i),  &
                npatch(1,1,i), npatch(1,2,i), npatch(1,3,i), npatch(1,4,i), &
                npatch(1,5,i), npatch(1,6,i), npatch(1,7,i),  &
                p23 &
             )

        ! 辺３の中点の曲面補間点
        eta=0.5; xi=0.5
        call fnpt_correct_pnt( &
                eta, xi,  &
                tri(1,1,i), tri(1,2,i), tri(1,3,i),  &
                npatch(1,1,i), npatch(1,2,i), npatch(1,3,i), npatch(1,4,i), &
                npatch(1,5,i), npatch(1,6,i), npatch(1,7,i),  &
                p31   &
             )

        ! ポリゴン１個目
        num_4 = num_4+1
        tri_4(1,1,num_4)=tri(1,1,i)
        tri_4(2,1,num_4)=tri(2,1,i)
        tri_4(3,1,num_4)=tri(3,1,i)
        tri_4(1,2,num_4)=p12(1)
        tri_4(2,2,num_4)=p12(2)
        tri_4(3,2,num_4)=p12(3)
        tri_4(1,3,num_4)=p31(1)
        tri_4(2,3,num_4)=p31(2)
        tri_4(3,3,num_4)=p31(3)

        ! ポリゴン２個目
        num_4 = num_4+1
        tri_4(1,1,num_4)=tri(1,2,i)
        tri_4(2,1,num_4)=tri(2,2,i)
        tri_4(3,1,num_4)=tri(3,2,i)
        tri_4(1,2,num_4)=p23(1)
        tri_4(2,2,num_4)=p23(2)
        tri_4(3,2,num_4)=p23(3)
        tri_4(1,3,num_4)=p12(1)
        tri_4(2,3,num_4)=p12(2)
        tri_4(3,3,num_4)=p12(3)

        ! ポリゴン３個目
        num_4 = num_4+1
        tri_4(1,1,num_4)=tri(1,3,i)
        tri_4(2,1,num_4)=tri(2,3,i)
        tri_4(3,1,num_4)=tri(3,3,i)
        tri_4(1,2,num_4)=p31(1)
        tri_4(2,2,num_4)=p31(2)
        tri_4(3,2,num_4)=p31(3)
        tri_4(1,3,num_4)=p23(1)
        tri_4(2,3,num_4)=p23(2)
        tri_4(3,3,num_4)=p23(3)

        ! ポリゴン４個目
        num_4 = num_4+1
        tri_4(1,1,num_4)=p12(1)
        tri_4(2,1,num_4)=p12(2)
        tri_4(3,1,num_4)=p12(3)
        tri_4(1,2,num_4)=p23(1)
        tri_4(2,2,num_4)=p23(2)
        tri_4(3,2,num_4)=p23(3)
        tri_4(1,3,num_4)=p31(1)
        tri_4(2,3,num_4)=p31(2)
        tri_4(3,3,num_4)=p31(3)
    enddo

    ! 面の法線ベクトル設定
    do i=1, num_4
        call get_tri_normal( tri_4(1,1,i), tri_4(1,2,i), tri_4(1,3,i), plane_norm4(1,i) )
    enddo

    ! 曲面補間した４倍のポリゴンSTLファイルに出力
    call output_stl_file( file_name_stl4_out, num_4, tri_4, plane_norm4 )

    write(*,*) '#### Program End'

end program
