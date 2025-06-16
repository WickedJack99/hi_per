#include <mpi.h>
#include "super_grid.h"

std::vector<MPI_Request> SuperGrid::receive_halos(HaloLayers& halo_layers)
{
    std::vector<MPI_Request> recv_requests(8);
    std::cout << halo_layers.top_halo.size() << "top halo size" << std::endl;
    MPI_Irecv(halo_layers.top_halo.data(), halo_layers.top_halo.size(), MPI_DOUBLE, neighbors_.top, neighbors_.top, comm_, &recv_requests[0]);
    MPI_Irecv(halo_layers.right_halo.data(), halo_layers.right_halo.size(), MPI_DOUBLE, neighbors_.right, neighbors_.right, comm_, &recv_requests[1]);
    MPI_Irecv(halo_layers.bottom_halo.data(), halo_layers.bottom_halo.size(), MPI_DOUBLE, neighbors_.bottom, neighbors_.bottom, comm_, &recv_requests[2]);
    MPI_Irecv(halo_layers.left_halo.data(), halo_layers.left_halo.size(), MPI_DOUBLE, neighbors_.left, neighbors_.left, comm_, &recv_requests[3]);

    MPI_Irecv(&halo_layers.top_right_corner, 1, MPI_DOUBLE, neighbors_.top_right, neighbors_.top_right, comm_, &recv_requests[4]);
    MPI_Irecv(&halo_layers.bottom_right_corner, 1, MPI_DOUBLE, neighbors_.bottom_right, neighbors_.bottom_right, comm_, &recv_requests[5]);
    MPI_Irecv(&halo_layers.top_left_corner, 1, MPI_DOUBLE, neighbors_.top_left, neighbors_.top_left, comm_, &recv_requests[6]);
    MPI_Irecv(&halo_layers.bottom_left_corner, 1, MPI_DOUBLE, neighbors_.bottom_left, neighbors_.bottom_left, comm_, &recv_requests[7]);

    return recv_requests;
}

std::vector<double> SuperGrid::get_inner_top_row()
{
    std::vector<double> row(cols());
    for (int j = 0; j < cols(); ++j)
        row[j] = (*this)(j, 0);
    return row;
}

std::vector<double> SuperGrid::get_inner_bottom_row()
{
    std::vector<double> row(cols());
    for (int j = 0; j < cols(); ++j)
        row[j] = (*this)(j, rows() - 1);
    return row;
}

std::vector<double> SuperGrid::get_inner_left_column()
{
    std::vector<double> col(rows());
    for (int i = 0; i < rows(); ++i)
        col[i] = (*this)(0, i);
    return col;
}

std::vector<double> SuperGrid::get_inner_right_column()
{
    std::vector<double> col(rows());
    for (int i = 0; i < rows(); ++i)
        col[i] = (*this)(cols() - 1, i);
    return col;
}

double SuperGrid::get_inner_top_left_corner()
{
    return (*this)(0, 0);
}

double SuperGrid::get_inner_top_right_corner()
{
    return (*this)(cols() - 1, 0);
}

double SuperGrid::get_inner_bottom_left_corner()
{
    return (*this)(0, rows() - 1);
}

double SuperGrid::get_inner_bottom_right_corner()
{
    return (*this)(cols() - 1, rows() - 1);
}

std::vector<MPI_Request> SuperGrid::inform_neighbors()
{
    std::vector<double> inner_top_row = get_inner_top_row();
    std::vector<double> inner_right_column = get_inner_right_column();
    std::vector<double> inner_bottom_row = get_inner_bottom_row();
    std::vector<double> inner_left_column = get_inner_left_column();

    double inner_top_right_corner = get_inner_top_right_corner();
    double inner_bottom_right_corner = get_inner_bottom_right_corner();
    double inner_top_left_corner = get_inner_top_left_corner();
    double inner_bottom_left_corner = get_inner_bottom_left_corner();

    if (rank_ == 0)
    {
        std::cout << "Rank " << rank_ << " neighbors:" << std::endl;
        std::cout << "  Top:          " << neighbors_.top << std::endl;
        std::cout << "  Right:        " << neighbors_.right << std::endl;
        std::cout << "  Bottom:       " << neighbors_.bottom << std::endl;
        std::cout << "  Left:         " << neighbors_.left << std::endl;
        std::cout << "  Top-Right:    " << neighbors_.top_right << std::endl;
        std::cout << "  Bottom-Right: " << neighbors_.bottom_right << std::endl;
        std::cout << "  Top-Left:     " << neighbors_.top_left << std::endl;
        std::cout << "  Bottom-Left:  " << neighbors_.bottom_left << std::endl;
    }

    std::vector<MPI_Request> send_requests(8);
    std::cout << inner_bottom_row.size() << "inner_bottom_row" << std::endl;
    MPI_Isend(inner_top_row.data(), inner_top_row.size(), MPI_DOUBLE, neighbors_.top, rank_, comm_, &send_requests[0]);
    MPI_Isend(inner_right_column.data(), inner_right_column.size(), MPI_DOUBLE, neighbors_.right, rank_, comm_, &send_requests[1]);
    MPI_Isend(inner_bottom_row.data(), inner_bottom_row.size(), MPI_DOUBLE, neighbors_.bottom, rank_, comm_, &send_requests[2]);
    MPI_Isend(inner_left_column.data(), inner_left_column.size(), MPI_DOUBLE, neighbors_.left, rank_, comm_, &send_requests[3]);

    MPI_Isend(&inner_top_right_corner, 1, MPI_DOUBLE, neighbors_.top_right, rank_, comm_, &send_requests[4]);
    MPI_Isend(&inner_bottom_right_corner, 1, MPI_DOUBLE, neighbors_.bottom_right, rank_, comm_, &send_requests[5]);
    MPI_Isend(&inner_top_left_corner, 1, MPI_DOUBLE, neighbors_.top_left, rank_, comm_, &send_requests[6]);
    MPI_Isend(&inner_bottom_left_corner, 1, MPI_DOUBLE, neighbors_.bottom_left, rank_, comm_, &send_requests[7]);

    return send_requests;
}

void SuperGrid::update()
{
    HaloLayers halo_layers(rows(), cols());
    std::vector<MPI_Request> send_requests = inform_neighbors();
    std::vector<MPI_Request> recv_requests = receive_halos(halo_layers);

    std::vector<MPI_Status> recv_status(recv_requests.size());
    MPI_Waitall(recv_requests.size(), recv_requests.data(), recv_status.data());

    std::vector<MPI_Status> send_status(send_requests.size());
    MPI_Waitall(send_requests.size(), send_requests.data(), send_status.data());

    merge_halos(halo_layers);
}

void SuperGrid::merge_halos(HaloLayers halo_layers)
{
    grid_(0, 0) = halo_layers.top_left_corner;
    grid_(grid_.cols() - 1, 0) = halo_layers.top_right_corner;
    grid_(0, grid_.rows() - 1) = halo_layers.bottom_left_corner;
    grid_(grid_.cols() - 1, grid_.rows() - 1) = halo_layers.bottom_right_corner;

    for (int col = 1; col < grid_.cols() - 1; col++)
    {
        grid_(col, 0) = halo_layers.top_halo[col - 1];
    }

    for (int col = 1; col < grid_.cols() - 1; col++)
    {
        grid_(col, grid_.rows() - 1) = halo_layers.bottom_halo[col - 1];
    }

    for (int row = 1; row < grid_.rows() - 1; row++)
    {
        grid_(0, row) = halo_layers.left_halo[row - 1];
    }

    for (int row = 1; row < grid_.rows() - 1; row++)
    {
        grid_(grid_.cols() - 1, row) = halo_layers.right_halo[row - 1];
    }
}

const Matrix SuperGrid::get_grid() const
{
    return grid_;
}